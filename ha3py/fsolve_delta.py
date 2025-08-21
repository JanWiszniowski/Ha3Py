r"""
The :math:`m_{max}` estimation by solving :math:`m_{max} = m_{max}^{obs} + \Delta`
--------------------------------------------------------------------------------------

:math:`m_{max}` is assessed by the solving the :math:`m_{max} = m_{max}^{obs} + \Delta` equation
by the *SciPy.fsolve* function.
The algorithm name in the configuration is 'solve_delta'

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
from ha3py.Delta import get_delta
from ha3py.utils import HaPyException


def _m_max_equation_to_solve(x, delta, pars):
    delta.m_max = x
    fun = x - pars['m_max_obs'] - delta(time=pars['time_span'], annual_lambda=pars['lambda_ref'])
    return fun


def m_max_solve_equation(configuration, magnitude_distribution=None, m_max=None, delta=None):
    r"""
    The function m_max_determination asses maximum regional magnitude by numerical solving
    the formula (Kijko, 1983, 1985; Pisarenko, 1991; Pisarenko et al., 1996)

    .. math::
        m_{max} - m_{max}^{obs}+\Delta = 0

        \sigma_{m_{max}} = \sqrt{\sigma_{m_{max}^{obs}}^2+\Delta^2}

    where :math`\Delta` is a class defined outside the function and the object ot the class is one
    of the Parameters of the call.

    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
        Required parameters for m_max_solve_equation in the dictionary are (keys are strings):

        * m_max_obs
        * sd_m_max_obs
        * m_min
        * m_max_current
        * time_span
        * theta coefficients, e.q. beta, lambda
        * constant coefficients, e,q, q_beta, q_lambda
        * m_max_current: starting m_max for solving the formula

    :type configuration: dict
    :param magnitude_distribution:
    :type magnitude_distribution:
    :param m_max:
    :type m_max:
    :param delta:
    :type delta:
    :return:
    :return:
        m_max       - estimated maximum magnitude
        sd_m_max    - standard deviation of maximum magnitude
    :rtype: (float, float)
    """
    if delta is None:
        delta = get_delta(configuration, magnitude_distribution=magnitude_distribution, m_max=m_max)
    m_max_obs = configuration['m_max_obs']
    sd_m_max_obs = configuration['sd_m_max_obs']
    root = fsolve(_m_max_equation_to_solve, np.array(m_max_obs), args=(delta, configuration))
    if not root:
        HaPyException(f'Solver can not find the solution for {delta.parameter_name}')
    if len(root) > 1:
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


if __name__ == "__main__":
    import sys
    import json

    if len(sys.argv) >= 2:
        with open(sys.argv[1], "r") as f:
            Configuration = json.load(f)
    else:
        Configuration = {'beta': 2.303, 'q_beta': 100.0, 'lambda': 2.0, 'q_lamb': 100.0, 'time_span': 500.0,
                         'm_min': 3.0, 'm_max_obs': 6.0, 'sd_m_max_obs': 0.25, 'm_max_current': 6.0,
                         'magnitude_distribution': 'Compound Gutenberg-Richter', 'delta': 'Kijko-Sellevoll'}
    # Using compound Gutenberg-Richter magnitude distribution and Kijko-Sellevoll method
    M_max, Sd_mag_max = m_max_solve_equation(Configuration)
    print(f'm_max = {M_max} (sd = {Sd_mag_max})')
