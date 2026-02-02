r"""
The :math:`m_{max}` estimation by the iteration
-----------------------------------------------------------------------------------------------

:math:`m_{max}` is assessed by the solving the :math:`m_{max} = m_{max}^{obs} + \Delta` equation
by the iteration.
The algorithm name in the configuration is 'iteration_delta'.

:copyright:
    Jan Wiszniowski <jwisz@igf.edu.pl>,
    Andrzej Kijko <andrzej.kijko@up.ac.za>
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
:version 0.0.1:
    2025-01-01
"""

from math import sqrt, fabs
from ha3py.Delta import get_delta


def m_max_solve_by_iteration(configuration, magnitude_distribution=None, m_max=None, delta=None):
    r"""
    The function m_max_determination asses maximum regional magnitude by the iteration
    according to the formula (Kijko, 1983, 1985; Pisarenko, 1991; Pisarenko et al., 1996)

    .. math::
        m_{max}^{(i+1)} = m_{max}^{obs}+\Delta_{m_{max}^{(i)}}

    where :math:`i` is the iteration step, :math:`m_{max}^{(0)}=m_{max}^{obs}`.
    Standard deviation is assumed as

    .. math::
        \sigma_{m_{max}} = \sqrt{\sigma_{m_{max}^{obs}}^2+\Delta^2}

    where the delta is a class defined outside the function and the object ot the class is one
    of the Parameters of the call.

    :param delta:
    :type delta:
    :param m_max:
    :type m_max:
    :param magnitude_distribution:
    :type magnitude_distribution:
    :param delta: Delta object (e.g. GuRiBaKiSe), which combines the magnitude distribution
        (e.g. Gutenberg-Richter-Bayes) and delta calculation method (e.g. Kijko-Sellovoll).
    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
        Required Parameters in the dictionary are (keys are strings):

            * m_max_obs
            * sd_m_max_obs
            * m_min
            * m_max_current
            * beta
            * lambda
            * time_span
            * q_beta
            * q_lamb
            * m_max_current: starting M_max in the iteration

    :type configuration: dict
    :return:
        M_max       - estimated maximum magnitude
        Sd_mag_max    - standard deviation of maximum magnitude

    """
    if delta is None:
        delta = get_delta(configuration, magnitude_distribution=magnitude_distribution, m_max=m_max)
    m_max_obs = configuration['m_max_obs']
    sd_m_max_obs = configuration['sd_m_max_obs']
    accuracy = 0.0001
    nr_iter_max = 20
    nr_iter = 0
    m_max_est_old = 0
    m_max_est_new = m_max_obs  # MAX OBS MAGNITUDE
    time = configuration['time_span']
    annual_lambda = configuration['lambda_ref']
    while (fabs(m_max_est_new - m_max_est_old) > accuracy) and (nr_iter_max > nr_iter):
        nr_iter += 1
        m_max_est_old = m_max_est_new
        delta.m_max = m_max_est_new
        m_max_est_new = m_max_obs + delta(time=time, annual_lambda=annual_lambda)
        print('Itr:{:d},M_max={}, delta={}'.format(nr_iter, m_max_est_new, type( delta).__name__))
    m_max = round(m_max_est_new, 2)
    if m_max <= m_max_obs:
        m_max = m_max_obs + .01
    if m_max > 9.99:
        m_max = 9.99
    sd_m_max = sqrt(sd_m_max_obs**2 + (m_max - m_max_obs)**2)
    sd_m_max = round(sd_m_max, 2)
    return m_max, sd_m_max


if __name__ == "__main__":
    import sys
    import json
    if len(sys.argv) >= 2:
        with open(sys.argv[1], "r") as f:
            parameters = json.load(f)
    else:
        parameters = {'beta': 2.303, 'q_beta': 100.0, 'lambda': 2.0, 'q_lamb': 100.0, 'time_span': 500.0,
                      'm_min': 3.0, 'm_max_obs': 6.0, 'sd_m_max_obs': 0.25, 'm_max_current': 6.0}
    # Using Gutenberg-Richter magnitude distribution and Kijko-Sellevoll delta method
    mag_max, sd_mag_max = m_max_solve_by_iteration(parameters)
    print('M_max = {} (sd = {})'.format(mag_max, sd_mag_max))
