r"""
Original Tate-Pisarenko :math:`m_{max}` assessment
---------------------------------------------------

The algorithm name in the configuration is 'Tate-Pisarenko'.
Do not miss with the Tate-Pisarenko by solve the equation.

The  :math:`m_{max}` is assessed by adding the :math:`\Delta=\frac{1}{nf\left( m_{max}^{obs} | m_{max}^{obs} \right)}`.

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

from ha3py.get_magnitude_distribution import get_magnitude_distribution
from ha3py.constant_values import EPS


def m_max_by_tate_pisarenko(configuration, magnitude_distribution=None, m_min=None):
    r"""
    The original Tate-Pisarenko method assumes the :math:`m_{max}` by formula.

    .. math::
        m_{max} = m_{max}^{obs} + \frac{1}{nf\left( m_{max}^{obs} | m_{max}^{obs} \right)}

        \sigma_{m_{max}} = \sigma_{m_{max}^{obs}} + \frac{n+1}{n^3 f^2 \left( m_{max}^{obs} | m_{max}^{obs} \right)}

    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict
    :param magnitude_distribution: Optional magnitude distribution object.
        If missing, the magnitude distribution object is created based on the configuration
    :type magnitude_distribution: MagnitudeDistribution
    :param m_min: Minimum value of the magnitude distribution.
        If missing, the maximum magnitude is taken from configuration
    :type m_min: float
    :return: Estimated maximum magnitude, standard deviation of maximum magnitude
    :rtype: (float, float)

    """
    m_max_obs = configuration['m_max_obs']
    sd_m_max_obs = configuration['sd_m_max_obs']
    n = configuration['time_span'] * configuration['lambda_ref']
    if magnitude_distribution is None:
        magnitude_distribution = get_magnitude_distribution(configuration, m_max=m_max_obs + EPS, m_min=m_min)
    else:
        magnitude_distribution = magnitude_distribution.copy()
        magnitude_distribution.m_max = m_max_obs + EPS
    f_m = magnitude_distribution(m_max_obs - EPS)
    return m_max_obs + 1.0 / (n * f_m), sd_m_max_obs + (n + 1) / (n ** 3 * f_m ** 2)
