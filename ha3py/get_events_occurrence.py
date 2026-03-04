from ha3py.PoissonOccurrence import PoissonOccurrence
from ha3py.GammaPoissonOccurrence import PoissonGammaCompoundOccurrence
from ha3py.utils import HaPyException


def get_events_occurrence(configuration, **kwargs):
    """
    Select and return the proper, defined in the configuration, events occurrence object.
    So far only Poisson and Poisson-gamma compound object can be defined.

    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict
    :param m_min: total minimum magnitude (optional, if it is missing, m_min taken from the configuration)
    :param m_max: Maximum value of the magnitude distribution.
        If missing, the maximum magnitude is taken from configuration
    :type m_max: float

    :return:
    """
    probability_name = configuration['occurrence_probability']
    if probability_name == 'Poisson':
        return PoissonOccurrence(configuration, **kwargs)
    elif probability_name == 'Poisson-gamma compound':
        return PoissonGammaCompoundOccurrence(configuration, **kwargs)
    else:
        raise HaPyException('Unknown events distribution')
