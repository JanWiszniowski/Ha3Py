"""
Ha3Py
(c) Jan Wiszniowski, Andrzej Kijko
ver. 2024-01
"""

import numpy as np
LN_10_ = 2.302585092994046
EPS = np.finfo(float).eps
EPS2 = 2.0 * EPS

METHODS = {
     1: 'Gibowicz-Kijko (1994)',
     2: 'Gibowicz-Kijko-Bayes (Kijko and Singh, 2011)',
     3: 'Kijko-Sellevoll (1989)',
     4: 'Kijko-Sellevoll-Bayes/compound (Kijko, 2004)',
     5: 'Tate-Pisarenko (Kijko and Graham, 1998)',
     6: 'Tate-Pisarenko-Bayes/compound (Kijko and Singh, 2011)',
     7: 'Non-Parametric-Gaussian/pseudo (Kijko, 2004)',
     8: 'Bayesian MEAN of shifted Likelihood Function & Gaussian Prior (Kijko, 2012)',
     9: 'Bayesian MEAN of Posterior Fiduicial & Prior Gauss (Kijko, 2004)',
     10:'Fixed maximum magnitude'
}

YMD_365 = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
YMD_366 = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
