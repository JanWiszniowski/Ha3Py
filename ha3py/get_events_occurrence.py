from ha3py.PoissonOccurrence import PoissonOccurrence
from ha3py.GammaPoissonOccurrence import PoissonGammaCompoundOccurrence
from ha3py.utils import HaPyException


def get_events_occurrence(configuration, **kwargs):
    probability_name = configuration['occurrence_probability']
    if probability_name == 'Poisson':
        return PoissonOccurrence(configuration, **kwargs)
    elif probability_name == 'Poisson-gamma compound':
        return PoissonGammaCompoundOccurrence(configuration, **kwargs)
    else:
        raise HaPyException('Unknown events distribution')
