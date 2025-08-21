"""
Module bayesian maximum magnitude computation
---------------------------------------------

..
    :copyright:
        Jan Wiszniowski <jwisz@igf.edu.pl>,
        Andrzej Kijko <andrzej.kijko@up.ac.za>
    :license:
        GNU Lesser General Public License, Version 3
        (https://www.gnu.org/copyleft/lesser.html)
    :version 0.0.1:
        2025-02-01

"""

import numpy as np
from ha3py.utils import HaPyException
import scipy.integrate as integrate

def bayesian_m_max(configuration, likelihood):
    r"""
    The bayesian_m_max function estimates the :math:`m_{max}` bze the Bayesian methods.
    The bayesian :math:`m_{max}` assessment methods incorporate any relevant information about :math:`m_{max}`
    which is defined as the apriori probability with the double-truncated Gaussian
    distribution :math:`\pi\left(m_{max}\right)=
    \mathcal{N}\left(m_{max}^{prior},\sigma_{m_{max}^{prior}}, m_{max}^L, m_{max}^U \right)`,
    where :math:`m_{\max}^U` is magnitude that might ever happen,
    and :math:`m_{max}^L` is usually the maximum observed magnitude :math:`m_{max}^{obs}`.
    Both the prior knowledge about :math:`m_{max}` and the knowledge coming from the observed magnitudes \mathbf{m},
    which is known as the posterior distribution of :math:`m_{max}`, summarise in the form

    .. math::

        p_{m_{max}}\left(m_{max}|\mathbf{m}\right)=
        \begin{cases}
        0 & : m_{max} < m_{\max}^L \\
        C\cdot \pi\left(m_{max}\right)\mathcal{L}\left(\mathbf{m}|m_{max}\right)
        & : m_{\max}^L \leqslant m_{max} \leqslant m_{\max}^U \\
        0 & : m_{max} > m_{\max}^U
        \end{cases},

    where :math:`C` is a normalising constant

    .. math::
        C=1/\int_{m_{max}^L}^{m_{max}^U}
        \pi\left(m_{max}\right)L\left(\mathbf{m}|m_{max}\right)dm_{max},

    and :math:`\mathcal{L}\left(\mathbf{m}|m_{max}\right)` is the likelihood function
    of measured earthquake magnitudes.
    Three Bayesian analogues of the maximum likelihood (ML) point estimators are used:
    the maximum posterior estimate (MAP) value, the posterior mean (PM) value (expected value),
    and the posterior median.

    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        including all catalogues
        and results of all computations.
    :type configuration: dict
    :param likelihood: The object describing the data depending on distribution
        :math:`\mathcal{L}\left(\mathbf{m}|m_{max}\right)`
        called the likelihood function,
    :type likelihood: BayesianBase
    :return: Estimated by bayesian method maximum magnitude
        and standard deviation of the maximum magnitude.
    :rtype: (float, float)
    """
    bayesian_m_max_estimator = configuration.get('bayesian_m_max_estimator', 'expected')
    if bayesian_m_max_estimator == 'expected':
        # return likelihood.expect()
        result = integrate.quad(lambda x: likelihood.pdf(x)*x, likelihood.m_max_l, likelihood.m_max_u)
        return result[0], likelihood.sd_m_max + likelihood.sd_prior_m_max
    elif bayesian_m_max_estimator == 'median':
        # return likelihood.median()
        space = np.linspace(likelihood.m_max_l, likelihood.m_max_u, 500)
        delta = (space[1]-space[0]) / 2.0
        aggregator = 0.0
        old_pdf = 0
        for idx, m in enumerate(space):
            new_pdf = likelihood.pdf(m)
            aggregator += (new_pdf + old_pdf) * delta
            old_pdf = new_pdf
            if aggregator >= 0.5:
                return m, likelihood.sd_m_max + likelihood.sd_prior_m_max
        raise HaPyException("Can not compute bayesian median")
    elif bayesian_m_max_estimator == 'max':
        max_pdf = 0
        m_max = None
        for m in np.linspace(likelihood.m_max_l, likelihood.m_max_u, 200):
            current_pdf = likelihood.pdf(m)
            if max_pdf < current_pdf:
                max_pdf = current_pdf
                m_max = m
        return m_max, likelihood.sd_m_max + likelihood.sd_prior_m_max
    else:
        raise HaPyException(f"Wrong bayesian m_max estimator ?'{bayesian_m_max_estimator}'")
