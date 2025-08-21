r"""
Non-bayesian Maximum possible earthquake magnitude estimation
-------------------------------------------------------------

:copyright:
    Jan Wiszniowski <jwisz@igf.edu.pl>,
    Andrzej Kijko <andrzej.kijko@up.ac.za>
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
:version 0.0.1:
    2025-01-01
"""

from ha3py.fsolve_delta import m_max_solve_equation
from ha3py.gibowicz_kijko import m_max_by_gibowicz_kijko
from ha3py.iteration import m_max_solve_by_iteration
from ha3py.momentum import m_max_by_momentum
from ha3py.primitive import m_max_primitive
from ha3py.robson_whitlock import m_max_by_robson_whitlock, m_max_by_robson_whitlock_cooke
from ha3py.utils import HaPyException


def non_bayesian_m_max_estimation(configuration, magnitude_distribution=None):
    r"""
    The estimation of :math:`\hat m_{max}` and :math:`\sigma_{m_{max}}`
    by the non-bayesian method defined in the configuration `m_max_assessment` parameter.
    The function chooses the method.

    :param magnitude_distribution:
    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict
    :return: Estimated by non bayesian method maximum magnitude
        and standard deviation of the maximum magnitude.
    :rtype: (float, float)
    """
    m_max_assessment = configuration.get('m_max_assessment', 'solve_delta')
    if m_max_assessment == 'solve_delta':
        return m_max_solve_equation(configuration, magnitude_distribution=magnitude_distribution)
    elif m_max_assessment == 'Gibowicz-Kijko':
        return m_max_by_gibowicz_kijko(configuration, magnitude_distribution=magnitude_distribution)
    elif m_max_assessment == 'iteration_delta':
        return m_max_solve_by_iteration(configuration, magnitude_distribution=magnitude_distribution)
    elif m_max_assessment == 'momentum':
        return m_max_by_momentum(configuration, magnitude_distribution=magnitude_distribution)
    elif m_max_assessment == 'primitive':
        return m_max_primitive(configuration)
    elif m_max_assessment == 'Robson-Whitlock':
        return m_max_by_robson_whitlock(configuration)
    elif m_max_assessment == 'Robson-Whitlock-Cooke':
        return m_max_by_robson_whitlock_cooke(configuration)
    else:
        raise HaPyException(f'Wrong m_max assessment procedure name ??{m_max_assessment}')
