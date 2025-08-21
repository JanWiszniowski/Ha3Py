r"""
Trivial :math:`m_{max}` assessment
----------------------------------

The algorithm name in the configuration is 'primitive'.

The primitive :math:`m_{max}` is assessed by adding the value 0.5 to :math:`m_{max}^{obs}`.

:copyright:
    Jan Wiszniowski <jwisz@igf.edu.pl>,
    Andrzej Kijko <andrzej.kijko@up.ac.za>
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
:version 0.0.1:
    2025-01-01
"""


def m_max_primitive(configuration):
    r"""
    The primitive method assumes the :math:`m_{max}` is greater than maximum observed magnitude by 0.5.

    .. math::
        m_{max} = m_{max}^{obs} + 0.5

        \sigma_{m_{max}} = \sigma_{m_{max}^{obs}}

    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict
    :return:
        m_max       - estimated maximum magnitude
        sd_m_max    - standard deviation of maximum magnitude
    :rtype: (float, float)

    """
    m_max_obs = configuration['m_max_obs']
    sd_m_max_obs = configuration['sd_m_max_obs']
    return m_max_obs + 0.5, sd_m_max_obs
