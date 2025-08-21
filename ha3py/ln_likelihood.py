"""
Likelihood coefficients estimation
----------------------------------

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
from ha3py.utils import HaPyException
from ha3py.get_events_occurrence import get_events_occurrence


def likelihood0(event_occurrence, catalogue):
    r"""
    Compute natural logarithm of likelihood for a catalogue with constant completeness level
    based on the non occurrence event of magnitude :math:`m` time distribution :math:`f_M^{max}`
    with the catalogue completeness magnitude

    .. math::
        ln\left( \mathcal{L}_\mathbf{\Theta} \right)=
        w_c\sum_{i=1}^{N}w_i\ln\left[f_M^{max} \left( m_i\middle|\mathbf{\Theta},t_i \right) \right],

    where :math:`N` is number of events in the catalogue, :math:`m_i` is the :math:`i`-th event magnitude,
    :param event_occurrence:
    :type event_occurrence:
    :math:`t_i` is the between event time, :math:`\mathbf{\Theta}` are the probability coefficients,
    :math:`w_c` is the weight of the current catalogue likelihood,
    and :math:`w_i` is the weight of the :math:`i`-th event

    :param catalogue: The catalogue must contain items:

        * 'm_min' (float) the completeness magnitude,
        * 'weight' (float, optional, default=1.0) the weight of the catalogue in the total ln likelihood,
            (:math:`w_c`)
        * 'events' (list) list of events. Each event is a dictionary containing:
            * 'magnitude' the event magnitude (:math:`m_i`),
            * 'time_span' time span the events (:math:`t_i`),
            * 'weight' (float, optional, default=1.0) the weight of the event in the catalogue ln likelihood,
                (:math:`w_i`).

    :type catalogue: dict
    :return: The likelihood of the catalogue
    :rtype: (float)

    Example, the gamma compound Poisson distribution case

    .. math::
        ln\left( \mathcal{L}_\lambda \mathcal{L}_\beta \right)=
        ln\left[ \lambda time\left( 1+\frac{\lambda time\left( 1-F_M\left( m \right) \right)}{q_\lambda} \right)
        ^{-\left( q_\lambda+1 \right)} \right]+ln\left( f_M\left( m \right) \right)

    """
    if not catalogue:
        return 0.0
    event_occurrence.m_min = catalogue['m_min']
    sum_ln_likelihood = 0.0
    for evt_phs in catalogue['earthquakes']:
        dtime = evt_phs['time_span']  # TIME INTERVAL [years]
        if dtime <= 0.0:
            dtime = 1.0
        magnitude = evt_phs['magnitude']  # MAGNITUDE OF PRE/HISTORIC EQ-s
        evt_weight = evt_phs.get('weight', 1.0)
        sum_ln_likelihood += event_occurrence.logpdf(magnitude, dtime) * evt_weight
    sum_ln_likelihood *= catalogue.get('weight', 1.0)
    return sum_ln_likelihood


def complete_catalogue_likelihood(event_occurrence, cat):
    if cat is None:
        return 0.0
    min_type = type(cat['m_min'])
    if min_type == float:
        return const_completeness_likelihood(event_occurrence, cat)
    elif min_type == str and type == 'time-varying':
        return time_varying_completeness_likelihood(event_occurrence, cat)
    else:
        raise HaPyException('Wrong m_min definition')


def const_completeness_likelihood(event_occurrence, catalogue):
    if not catalogue:
        return 0.0
    event_occurrence.m_min = catalogue['m_min']
    dtime = catalogue['time_span']
    n = len(catalogue['earthquakes'])
    if dtime <= 0.0:
        dtime = 1.0
    sum_ln_likelihood = event_occurrence.ln_d_pmf(n, dtime)
    # sum_ln_likelihood = np.log(event_occurrence.lamb) * n - np.log(dtime * event_occurrence.lamb + event_occurrence.q_lambda) * (event_occurrence.q_lambda + n)
    for evt_phs in catalogue['earthquakes']:
        evt_weight = evt_phs.get('weight', 1.0)
        magnitude = evt_phs['magnitude']
        sum_ln_likelihood += event_occurrence.magnitude_distribution.logpdf(magnitude) * evt_weight
    sum_ln_likelihood *= catalogue.get('weight', 1.0)
    return sum_ln_likelihood


def time_varying_completeness_likelihood(event_occurrence, catalogue):
    if not catalogue:
        return 0.0
    sum_ln_likelihood = 0.0
    for evt_phs in catalogue['earthquakes']:
        event_occurrence.m_min = evt_phs['m_min']
        dtime = evt_phs['time_span']
        if dtime <= 0.0:
            dtime = 1.0
        magnitude = evt_phs['magnitude']
        evt_weight = evt_phs.get('weight', 1.0)
        sum_ln_likelihood += event_occurrence.logpdf(magnitude, dtime)
        # sum_ln_likelihood += event_occurrence.magnitude_distribution.logpdf(magnitude) ????
        sum_ln_likelihood *= evt_weight  # LN[PDF(MAG)]
    sum_ln_likelihood *= catalogue.get('weight', 1.0)
    return sum_ln_likelihood


def ln_likelihood(x_v, configuration, m_max=None):
    """

    :param x_v:
        Array of event distribution and occurrence coefficients.
        They must agree with the definition of distributions
    :type x_v: numpy.array
    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict
    :param m_max: Maximum magnitude. If missing the maximum magnitude is taken from configuration
    :type m_max: float
    :return: The natural logarithm of likelihood computed ever all catalogues.
    :rtype: (float)

    """
    event_occurrence = get_events_occurrence(configuration, theta=x_v, m_max=m_max)
    part_phs = likelihood0(event_occurrence, configuration.get('paleo_catalog'))
    part_his = likelihood0(event_occurrence, configuration.get('historic_catalog'))
    part_const = 0.0  # CONSTRAINS
    part_prior = 0.0
    part_comp = 0.0
    p_comp = configuration.get('complete_catalogs')
    if p_comp:
        for complete_catalog in p_comp:
            part_comp += complete_catalogue_likelihood(event_occurrence, complete_catalog)
    coefficient_names = event_occurrence.coefficient_names[:-1]
    coefficient_values = event_occurrence.coefficients[:-1]
    sd_const = 0.1
    n_const = 4.0
    for idx, coefficient_name in enumerate(coefficient_names):
        coefficient_value = coefficient_values[idx]
        part_prior = - 0.5 * abs(coefficient_value - configuration.get(f'prior_{coefficient_name}',
                                                                       coefficient_value))
        part_prior /= (configuration.get(f'sd_prior_{coefficient_name}', 1.0)) ** 2
        min_value = configuration.get(f'min_{coefficient_name}')
        max_value = configuration.get(f'max_{coefficient_name}')
        if min_value is not None and coefficient_value < min_value:
            part_const = - ((min_value - coefficient_value) / sd_const) ** n_const
        elif max_value is not None and coefficient_value > max_value:
            part_const = - ((coefficient_value - max_value) / sd_const) ** n_const
    retval = - (part_phs + part_his + part_comp + part_prior + part_const)
    return retval
