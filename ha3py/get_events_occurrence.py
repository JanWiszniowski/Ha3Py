from ha3py.PoissonOccurrence import PoissonOccurrence
from ha3py.GammaPoissonOccurrence import PoissonGammaCompoundOccurrence
from ha3py.utils import HaPyException


def get_events_occurrence(configuration, **kwargs):
    """
    Select and return the proper, defined in the configuration, events occurrence object.
    So far only Poisson and Poisson-gamma compound object can be defined.

    :param configuration:
    :param m_min: total minimum magnitude (optional, if it is missing, m_min taken from the configuration)
    :param m_max: maximum magnitude
        (optional, if it is missing, proper m_max or current_m_max is taken from the configuration)
    :param theta: list of events occurrence coefficients - first lamba coefficients, next beta coefficients
        (optional, if it is missing, coefficients, like lambda, beta, are taken from the configuration)

    :return:
    """
    probability_name = configuration['occurrence_probability']
    if probability_name == 'Poisson':
        return PoissonOccurrence(configuration, **kwargs)
    elif probability_name == 'Poisson-gamma compound':
        return PoissonGammaCompoundOccurrence(configuration, **kwargs)
    else:
        raise HaPyException('Unknown events distribution')
