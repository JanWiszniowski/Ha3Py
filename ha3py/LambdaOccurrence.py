"""

..
    :copyright:
        Jan Wiszniowski <jwisz@igf.edu.pl>,
        Andrzej Kijko <andrzej.kijko@up.ac.za>
    :license:
        GNU Lesser General Public License, Version 3
        (https://www.gnu.org/copyleft/lesser.html)
    :version 0.0.1:
        2025-01-01


"""

from abc import ABC
from ha3py.BaseOccurrence import OccurrenceBase


class LambdaOccurrence(OccurrenceBase, ABC):
    """
    Base class of event occurrence classes, which are described by lambda coefficient.

    """
    def __init__(self, configuration, name, **kwargs):
        theta = kwargs.get('theta')
        if theta is None:
            self.lamb = configuration.get('lambda', 1.0)
        else:
            self.lamb = theta[0]
        self.lambda_origin = self.lamb
        self.m_min_origin = 1.0
        super().__init__(configuration, name, **kwargs)
        self.m_min_origin = self.magnitude_distribution.m_min

    def _set_m_min(self, val):
        self.magnitude_distribution.m_min = self.m_min_origin
        self.lamb = self.lambda_origin * self.magnitude_distribution.sf(val)
        # self.a = self.magnitude_distribution.m_min

    def _coefficient_names(self):
        return ['lambda']

    def _coefficient_values(self):
        return [self.lamb]

    def _d_mean(self, t):
        return self.lamb * t
