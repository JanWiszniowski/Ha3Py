r"""
The :math:`m_{max}` estimation by momentum
------------------------------------------

The algorithm name in the configuration is 'momentum'.

:copyright:
    Jan Wiszniowski <jwisz@igf.edu.pl>,
    Andrzej Kijko <andrzej.kijko@up.ac.za>
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
:version 0.0.1:
    2025-01-01
"""


import numpy as np
from math import sqrt
from ha3py.utils import HaPyException


def m_max_by_momentum_compute(mag_a):
    """

    :param mag_a:
    :return:
    """
    m_shift = mag_a - min(mag_a)  # work with distribution shifted to start support at zero
    n_mag = mag_a.size
    mom1 = sum(m_shift) / n_mag  # first moment
    mom2 = sum(m_shift ** 2) / n_mag  # second moment
    mom3 = sum(m_shift ** 3) / n_mag  # third moment
    a1 = mom2 - 2 * mom1 ** 2
    a2 = 3 * mom1 * mom2 - mom3
    a3 = 2 * mom1 * mom3 - 3 * mom2 ** 2
    d = sqrt(a2 ** 2 - 4 * a1 * a3)  # discriminant in estimator <m_max> as
    answer1 = (-a2 + d) / 2 / a1
    answer2 = (-a2 - d) / 2 / a1
    return max(answer1, answer2) + min(mag_a)


def sd_m_max_by_momentum_compute(mag_v):
    r"""
    function var_m_max = f_boot_var(mag_v)

    Subfunction f_boot_var calculates variance of moment estimator of M_max by
    means of bootstrap resampling

    :param mag_v: Vector of seismic event magnitudes
    :return: Estimated variance from bootstrap procedure

    """
    cat_size = mag_v.size
    boot_val = []
    sample_size = min(100, np.ceil(cat_size / 2))
    # boot_size = max(np.ceil(cat_size * 0.75), cat_size-10)
    for _ in range(sample_size):
        boot_sample = np.random.choice(mag_v, size=30, replace=False)
        boot_val.append(m_max_by_momentum_compute(boot_sample))
    return np.std(boot_val)


def get_magnitudes_for_momentum(configuration, magnitude_distribution=None):
    p_phs = configuration.get('paleo_catalog')
    p_his = configuration.get('historic_catalog')
    p_comp = configuration.get('complete_catalogs')
    if p_phs or p_his or not p_comp:
        HaPyException('Moment estimator can be applied only to complete_catalogs')
        return None
    mag_mom = []
    if magnitude_distribution is None:
        m_min = configuration['m_min']
    else:
        m_min = magnitude_distribution.m_min
    for complete_catalog in p_comp:
        mag_mom.extend([eq['magnitude'] for eq in complete_catalog['earthquakes'] if eq['magnitude'] >= m_min])
    return np.ndarray(mag_mom)


def m_max_by_momentum(configuration, magnitude_distribution=None):  # magnitude_distribution=None, m_max=None, delta=None
    r"""
     :math:`m_{max}` evaluation according to moment estimator (SEE "Estimation of
     Paramers of a Right Truncated Exponential Distribution" by U.J. Dixit and P.N. Nasiri
     published in Statistical Papers Vol 49 (2008) pp.225-236)

    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict
    :param magnitude_distribution: object
    :type magnitude_distribution: MagnitudeDistribution
    :return:
    :rtype:

    Function history:

        * Created by P.J. Vermeulen on FEB 2014
            (Created as additional procedure for the program mmax.m)
        * MAR 2014: Calculation of standard deviation by bootstrap method
            included. Done by means of subfunction f_boot_var.
        * 2025.01.01 code in Python

    """
    mag_a = get_magnitudes_for_momentum(configuration, magnitude_distribution=magnitude_distribution)
    if not mag_a:
        return None, None
    m_max = m_max_by_momentum_compute(mag_a)
    var_m_max = sd_m_max_by_momentum_compute(mag_a)
    return m_max, var_m_max
