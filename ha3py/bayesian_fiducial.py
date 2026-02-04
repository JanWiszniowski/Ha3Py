"""
Module bayesian fiducial maximum magnitude computation
------------------------------------------------------

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
from ha3py.m_max_utils import non_bayesian_m_max_estimation
from ha3py.BayesianMmax import BayesianBase
from ha3py.get_events_occurrence import get_events_occurrence


class BayesianFiducial(BayesianBase):
    r"""
    The class that estimate the :math:`m_{max}` by the Bayesian fiducial method.
    In the case of single catalog, the class define the likelihood function


    .. math::
        \mathcal{L}\left(\mathbf{m}|m_{max}\right)=
        f_{M_{max}}^{FID} \left( m_{max} \right)

    used in the Bayesian posterior distribution of :math:`m_{max}`

    .. math::
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
        \pi\left(m_{max}\right)L\left(\mathbf{m}|m_{max}\right)dm_{max},

    and :math:`m_{max}^L=m_{max}^{obs}`.

    The probability density function is

    .. math::
        f_{M_{max}}^{FID} \left( m_{max} \right)=
        n\left[F_M\left(m_{max}^{obs}|m_{max} \right)\right]^{n-1}
        \frac{\partial S_M\left(m_{max}^{obs}|m_{max} \right)}{\partial m_{max}}

    In the case of single catalog cumulate density function is

    .. math::
        F_{M_{max}}^{FID}(m_{max})=
        1-\left[F_M\left(m_{max}^{obs}|m_{max} \right)\right]^n

    In the case of many catalogs

    .. math::
        F_{M_{max}}^{FID}(m_{max})=
        1-\prod_{i}^{k}\left[F_M\left(m_{max}^{obs}|m_{max},m_{min}^{(i)} \right)\right]^{n_i},

    where :math:'k' is the number of catalogues,
    :math:'m_{min}^{(i)}' is the completeness level in the catalogue,
    and :math:'n_i' is the number of events in the catalogue.

    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict


    """
    def __init__(self, configuration, magnitude_distribution=None):
        """

        :param configuration:
        :type configuration:
        """
        m_max, sd_m_max = non_bayesian_m_max_estimation(configuration)
        fiducial_m_min = configuration.get('fiducial_m_min')
        if fiducial_m_min is None:
            event_occurrence = get_events_occurrence(configuration, m_max=m_max)
        else:
            event_occurrence = get_events_occurrence(configuration, m_max=m_max, m_min=fiducial_m_min)
        self.n = event_occurrence.d_expected(t=configuration['time_span'])
        self.magnitude_distribution = event_occurrence.magnitude_distribution
        self.m_max_obs = configuration['m_max_obs']
        super().__init__(configuration, magnitude_distribution=magnitude_distribution,
                         m_max_pair=(m_max, sd_m_max))

    def _likelihood(self, m_max):
        """

        :param m_max:
        :type m_max:
        :return:
        :rtype:
        """
        if np.isscalar(m_max):
            self.magnitude_distribution.m_max = m_max
            o = float(self.n)
            o *= self.magnitude_distribution.cdf(self.m_max_obs) ** (self.n - 1)
            o *= self.magnitude_distribution.grad_sf(self.m_max_obs, coefficient_name='m_max')
            return o
        else:
            n = m_max.size
            output = np.zeros(n)
            for idx in range(n):
                self.magnitude_distribution.m_max = m_max[idx]
                o = float(self.n)
                o *= self.magnitude_distribution.cdf(self.m_max_obs) ** (self.n - 1)
                o *= self.magnitude_distribution.grad_sf(self.m_max_obs, coefficient_name='m_max')
                output[idx] = o
            return output


def get_bayesian_fiducial(configuration, magnitude_distribution=None):
    """

    :param configuration:
    :type configuration:
    :param magnitude_distribution:
    :type magnitude_distribution: object
    :return:
    :rtype:
    """
    return BayesianFiducial(configuration, magnitude_distribution=magnitude_distribution)
