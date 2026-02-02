r"""
Module bayesian maximum magnitude computation assuming normal distribution
--------------------------------------------------------------------------

In this method we assume the normal distribution of a prior maximum magnitude
:math:`\mathcal{N}\left( m_{max}^{prior}, \sigma_{m_{max}^{prior}} \right)`
and normal distribution of a maximum magnitude estimated based on catalogues
:math:`\mathcal{N}\left( m_{max}, \sigma_{m_{max}} \right)`.


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
from scipy.stats import truncnorm
from ha3py.m_max_utils import non_bayesian_m_max_estimation
from ha3py.BayesianMmax import init_bayesian_m_max


def m_max_bayesian_norm(configuration, magnitude_distribution=None):
    r"""
    The m_max_by_bayesian_norm assumes the normal distribution of a prior maximum magnitude
    :math:`\mathcal{N}\left( m_{max}^{prior}, \sigma_{m_{max}^{prior}} \right)`
    and normal distribution of a maximum magnitude :math:`\mathcal{N}\left( m_{max}, \sigma_{m_{max}} \right)`
    estimated based on catalogues.
    Coefficients :math:`m_{max}` and :math:`\sigma_{m_{max}}` are assesses by one of non bayesian method.
    The posterior maximum magnitude and its standard deviation are then

    .. math::

        m^{posterior}_{max}=
        \frac{{\sigma^2_{m_{max}^{prior}}m_{max}}+\sigma^2_{m_{max}}m_{\max}^{prior}}
        {\sigma^2_{m_{max}^{prior}}+\sigma^2_{m_{max}}},

        \sigma_{m_{max}^{posterior}}=
        \frac{\sigma_{m_{max}^{prior}}\sigma_{m_{max}}}
        {\sqrt{\sigma_{m_{max}^{prior}}+\sigma_{m_{max}}}}

    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict
    :param magnitude_distribution:
    :type magnitude_distribution: object
    :return: The posterior maximum magnitude and its standard deviation.
    :rtype: tuple(float, float)
    """
    # m_max, sd_m_max = non_bayesian_m_max_estimation(configuration, magnitude_distribution=magnitude_distribution)
    # # m_max, sd_m_max = 1, 2
    # prior_m_max = configuration['prior_m_max']
    # sd_prior_m_max = configuration['sd_prior_m_max']
    m_max, sd_m_max, prior_m_max, sd_prior_m_max = init_bayesian_m_max(
        configuration, magnitude_distribution=magnitude_distribution)
    d = sd_m_max ** 2 + sd_prior_m_max ** 2
    return ((sd_m_max ** 2 * prior_m_max + sd_prior_m_max ** 2 * m_max) / d,
            sd_m_max * sd_prior_m_max / np.sqrt(d))


def get_bayesian_truncnorm(configuration, magnitude_distribution=None):
    r"""
    The m_max_by_bayesian_norm assumes the truncated normal distribution of a prior maximum magnitude
    :math:`\mathcal{N}\left( m_{max}^{prior},
    \sigma_{m_{max}^{prior}}, m_{max}^{L}, m_{max}^{U} \right)`
    and normal distribution of a maximum magnitude
    :math:`\mathcal{N}\left( m_{max}, \sigma_{m_{max}} \right)`
    estimated based on catalogues,
    where :math:`m_{max}^{U}` is the maximum possible magnitude that we believe might ever happen and
    :math:`m_{max}^{L}` is the magnitude that we are sure the maximum magnitude is greater.
    Coefficients :math:`m_{max}` and :math:`\sigma_{m_{max}}` are assesses by one of non bayesian method.

    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict
    :param magnitude_distribution:
    :type magnitude_distribution:
    :return: The
    :rtype: tuple(float, float)
    """
    m_max, sd_m_max, prior_m_max, sd_prior_m_max = init_bayesian_m_max(
        configuration, magnitude_distribution=magnitude_distribution)
    m_max_u = configuration.get('upper_m_max', 9.5)
    m_max_obs = configuration['m_max_obs']
    m_max_l = configuration.get('lower_m_max', m_max_obs)
    d = sd_m_max ** 2 + sd_prior_m_max ** 2
    loc = (sd_m_max ** 2 * prior_m_max + sd_prior_m_max ** 2 * m_max) / d
    scale = sd_m_max * sd_prior_m_max / np.sqrt(d)
    a = (m_max_l - loc) / scale
    b = (m_max_u - loc) / scale
    distribution = truncnorm(a, b, loc=loc, scale=scale)
    distribution.m_max_l = m_max_l
    distribution.m_max_u = m_max_u
    distribution.sd_m_max = sd_m_max
    distribution.sd_prior_m_max = sd_prior_m_max
    return distribution
    # m, v, _, _ = truncnorm.stats(m_max_l, m_max_u,
    #                              loc=(sd_m_max ** 2 * prior_m_max + sd_prior_m_max ** 2 * m_max) / d,
    #                              scale=sd_m_max * sd_prior_m_max / np.sqrt(d))
    # return m, np.sqrt(v)
