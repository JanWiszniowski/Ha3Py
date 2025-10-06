r"""
Maximum possible earthquake magnitude estimation
------------------------------------------------

:copyright:
    Jan Wiszniowski <jwisz@igf.edu.pl>,
    Andrzej Kijko <andrzej.kijko@up.ac.za>
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
:version 0.0.1:
    2025-01-01
"""

from ha3py.bayesian_estimators import bayesian_m_max
from ha3py.bayesian_fiducial import get_bayesian_fiducial
from ha3py.bayesian_by_shift import get_bayesian_by_shift
from ha3py.m_max_utils import non_bayesian_m_max_estimation
from ha3py.bayesian_normal import m_max_bayesian_norm, get_bayesian_truncnorm
from ha3py.get_events_occurrence import get_events_occurrence
from ha3py.utils import HaPyException
from ha3py.configuration import load_configuration, save_configuration


def m_max_estimation(configuration):
    r"""
    The estimation of :math:`\hat m_{max}` and :math:`\sigma_{m_{max}}`
    by the one of all available methods, both bayesian and non-bayesian.
    The non-bayesian method is defined in the configuration "m_max_assessment" parameter.
    The bayesian method is defined in the configuration "bayesian_m_max_assessment" parameter.
    If the "bayesian_m_max_assessment" parameter is missing or empty,
    the :math:`\hat m_{max}` is estimated by the non-bayesian method.

    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict
    :return: Estimated by non bayesian method maximum magnitude
        and standard deviation of the maximum magnitude.
    :rtype: (float, float)

    """
    event_occurrence = get_events_occurrence(configuration)
    if 'm_min_ref' in configuration:
        event_occurrence.m_min = configuration['m_min_ref']
    configuration['lambda_ref'] = event_occurrence.d_mean()

    bayesian_assessment = configuration.get('bayesian_m_max_assessment')
    if not bayesian_assessment:
        return non_bayesian_m_max_estimation(configuration, magnitude_distribution=event_occurrence.magnitude_distribution)
    elif bayesian_assessment == 'bayesian_normal':
        return bayesian_m_max(configuration, get_bayesian_truncnorm(configuration, magnitude_distribution=event_occurrence.magnitude_distribution))
    elif bayesian_assessment == 'bayesian_normal_unlimited':
        return m_max_bayesian_norm(configuration, magnitude_distribution=event_occurrence.magnitude_distribution)
    elif bayesian_assessment == 'bayesian_by_shift':
        return bayesian_m_max(configuration, get_bayesian_by_shift(configuration, magnitude_distribution=event_occurrence.magnitude_distribution))
    elif bayesian_assessment == 'bayesian_fiducial':
        return bayesian_m_max(configuration, get_bayesian_fiducial(configuration, magnitude_distribution=event_occurrence.magnitude_distribution))
    elif bayesian_assessment == 'fixed value':
        return configuration['prior_m_max'], configuration['sd_prior_m_max']
    else:
        raise HaPyException(f"Wrong bayesian m_max assessment procedure name ??{bayesian_assessment}")


def main():
    params = load_configuration()
    m_max, sd_m_max = m_max_estimation(params)
    print(f"Result of m_max estimation: m_max = {m_max:4.2f} (+/- {sd_m_max:4.2f})")
    if input('Save the result? [yes/no] >') == 'yes':
        params['m_max'] = m_max
        params['sd_m_max'] = sd_m_max
        save_configuration(params)


if __name__ == "__main__":
    main()
