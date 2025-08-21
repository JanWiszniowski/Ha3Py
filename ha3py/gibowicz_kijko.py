r"""
The :math:`m_{max}` assessment by Gibowicz-Kijko procedure
----------------------------------------------------------

The algorithm name in the configuration is 'Gibowicz-Kijko'.

:copyright:
    Jan Wiszniowski <jwisz@igf.edu.pl>,
    Andrzej Kijko <andrzej.kijko@up.ac.za>
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
:version 0.0.1:
    2025-01-01

"""

from math import sqrt
from scipy.optimize import fsolve
import numpy as np
from ha3py.utils import HaPyException
from ha3py.get_magnitude_distribution import get_magnitude_distribution


def _m_max_equation_to_solve(x, magnitude_distribution, n, m_max_obs):
    magnitude_distribution.m_max = x
    return magnitude_distribution.cdf(m_max_obs) - n / (n+1)


def m_max_by_gibowicz_kijko(configuration, magnitude_distribution=None, m_max=None, m_min=None):
    r"""
    Gibowicz-Kijko procedure estimates :math:`m_{max}` by numerically solving the equation

    .. math::
        F_M\left( m_{max}^{obs} \right) - \frac{n}{n+1} = 0,

    where :math:`n` can be assumed as

    .. math::
        n=\lambda time.

    Standard deviation is assumed as

    .. math::
        \sigma_{m_{max}} = \sqrt{\sigma_{m_{max}^{obs}}^2+(m_{max}-m_{max}^{obs})^2}

    :param m_min:
    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict
    :param m_max:
    :type m_max:
    :param magnitude_distribution:
    :type magnitude_distribution:
    :return: estimated maximum magnitude, standard deviation of maximum magnitude.
    :rtype: (float, float)

    """
    if magnitude_distribution is None:
        magnitude_distribution = get_magnitude_distribution(configuration, m_max=m_max, m_min=m_min)
    n = configuration['time_span'] * configuration['lambda_ref']
    m_max_obs = configuration['m_max_obs']
    sd_m_max_obs = configuration['sd_m_max_obs']
    root = fsolve(_m_max_equation_to_solve, np.array(m_max_obs),
                  args=(magnitude_distribution, n, m_max_obs))
    if not root:
        HaPyException(f'G-K Solver can not find the solution for {magnitude_distribution.parameter_name}')
    if len(root) > 1:
        print(f'WARNING! G-KSolver found for {magnitude_distribution.parameter_name} many solutions {root}')
        m_max = np.max(root)
    else:
        m_max = root[0]
    if m_max <= m_max_obs:
        m_max = m_max_obs + 0.01
    if m_max > 9.99:
        m_max = 9.99
    sd_m_max = sqrt(sd_m_max_obs ** 2 + (m_max - m_max_obs) ** 2)
    sd_m_max = round(sd_m_max, 2)
    return m_max, sd_m_max
