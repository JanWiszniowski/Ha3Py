r"""
:math:`m_{max}` assessment by Robson-Whitlock and Robson-Whitlock-Cooke procedure
---------------------------------------------------------------------------------

:math:`m_{max}` is assessed by adding to :math:`m_{max}^{obs}` the difference
between :math:`m_{max}^{obs}` and the second maximum magnitude.

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


def get_magnitudes(configuration):
    m_n = -100.0
    m_n_1 = None
    p_phs = configuration['paleo_catalog']
    p_his = configuration['historic_catalog']
    p_comp = configuration['complete_catalogs']
    if p_phs:
        for eq in p_phs['earthquakes']:
            if eq['magnitude'] > m_n:
                m_n_1 = m_n
                m_n = eq['magnitude']
    if p_his:
        for eq in p_his['earthquakes']:
            if eq['magnitude'] > m_n:
                m_n_1 = m_n
                m_n = eq['magnitude']
    if p_comp:
            for complete_catalog in p_comp:
                for eq in complete_catalog:
                    if eq['magnitude'] > m_n:
                        m_n_1 = m_n
                        m_n = eq['magnitude']
    return m_n_1


def m_max_by_robson_whitlock(configuration):
    r"""
    Robson-Whitlock procedure

    .. math::
        m_{max} = m_{max}^{obs}+\left( m_{max}^{obs} - m_{n-1} \right)

        \sigma_{m_{max}} = \sqrt{5\sigma_{m_{max}^{obs}}^2 + \left( m_{max}^{obs} - m_{n-1} \right)^2}

    The algorithm name in the configuration is 'Robson-Whitlock'.

    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict
    :return:
        m_max       - estimated maximum magnitude
        sd_m_max    - standard deviation of maximum magnitude
    :rtype: (float, float)

    """
    m_n_1 = get_magnitudes(configuration)
    m_max_obs = configuration['m_max_obs']
    m_max = 2.0 * m_max_obs - m_n_1
    sd_m_max_obs = configuration['sd_m_max_obs']
    sd_m_max = sqrt(5.0 * sd_m_max_obs ** 2 + (m_max_obs - m_n_1) ** 2)
    return m_max, sd_m_max


def m_max_by_robson_whitlock_cooke(configuration):
    r"""
    Robson-Whitlock-Cooke procedure

    .. math::
        m_{max} = m_{max}^{obs}+0.5\left( m_{max}^{obs} - m_{n-1} \right)

        \sigma_{m_{max}} = \sqrt{1.5\sigma_{m_{max}^{obs}}^2 + 0.25\left( m_{max}^{obs} - m_{n-1} \right)^2}

    The algorithm name in the configuration is 'Robson-Whitlock-Cooke'.


    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict
    :return:
        m_max       - estimated maximum magnitude
        sd_m_max    - standard deviation of maximum magnitude
    :rtype: (float, float)

    """
    m_n_1 = get_magnitudes(configuration)
    m_max_obs = configuration['m_max_obs']
    m_max = m_max_obs + 0.5 * (m_max_obs - m_n_1)
    sd_m_max_obs = configuration['sd_m_max_obs']
    sd_m_max = sqrt(1.5 * sd_m_max_obs ** 2 + 0.25 * (m_max_obs - m_n_1) ** 2)
    return m_max, sd_m_max
