r"""
The :math:`m_{max}` assessment by direct solution of Tate-Pisarenko
-------------------------------------------------------------------

The solution requires the Gutenberg-Richter magnitude distribution

..
    :copyright:
        Jan Wiszniowski <jwisz@igf.edu.pl>,
        Andrzej Kijko <andrzej.kijko@up.ac.za>
    :license:
        GNU Lesser General Public License, Version 3
        (https://www.gnu.org/copyleft/lesser.html)
    :version 0.0.1:
        2026-08-03

"""

from math import sqrt
from scipy.optimize import fsolve
from scipy import lambertw
import numpy as np
from ha3py.utils import HaPyException
from ha3py.get_magnitude_distribution import get_magnitude_distribution



def m_max_by_direct_tate_pisarenko(configuration, magnitude_distribution=None, m_max=None, m_min=None):
    r"""
    :math:`m_{max}`, which is the direct solution of Tate-Pisarenko, is described by the solution

    .. math::
        {\hat{m}}_{max}=m_{max}^{obs}-
        \frac{W\left[-\frac{exp\left(\beta m_{max}^{obs}-\beta m_{min}-\frac{1}{n}\right)}{n}\right]}{\beta}
        -\frac{1}{n\beta}

    where :math:`W` is the Lambert W function and :math:`n=\lambda time`.

    Standard deviation is assumed as

    .. math::
        \sigma_{m_{max}} = \sqrt{\sigma_{m_{max}^{obs}}^2+(m_{max}-m_{max}^{obs})^2}

    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict
    :param m_max: Maximum value of the magnitude distribution.
        If missing, the maximum magnitude is taken from configuration
    :type m_max: float
    :param m_min: Minimum value of the magnitude distribution.
        If missing, the maximum magnitude is taken from configuration
    :type m_min: float
    :param magnitude_distribution:  Optional magnitude distribution object.
        If missing, the magnitude distribution object is created based on the configuration
    :type magnitude_distribution: MagnitudeDistribution
    :return: estimated maximum magnitude, standard deviation of maximum magnitude.
    :rtype: (float, float)

    """
    distribution_name = configuration.get('magnitude_distribution', 'Gutenberg-Richter')
    if distribution_name != 'Gutenberg-Richter':
        print('The direct Tate-Pisarenko requires Gutenberg-Richter magnitude distribution')
        return None, None
    n = configuration['time_span'] * configuration['lambda_ref']
    beta = configuration['beta']
    m_max_obs = configuration['m_max_obs']
    sd_m_max_obs = configuration['sd_m_max_obs']
    if m_max_obs - m_min - np.log(n) / beta + (n - 1) / n / beta > 0:
        print(f"The direct Tate-Pisarenko requires n >= {np.exp(beta * (m_max_obs- m_min) + 1)}")
        return None, None
    m_max = m_max_obs - lambertw(-np.exp(beta * (m_max_obs - m_min) + 1 / n) / n) / beta - 1 / n / beta
    if m_max <= m_max_obs:
        m_max = m_max_obs + 0.01
    if m_max > 9.99:
        m_max = 9.99
    sd_m_max = sqrt(sd_m_max_obs ** 2 + (m_max - m_max_obs) ** 2)
    sd_m_max = round(sd_m_max, 2)
    return m_max, sd_m_max
