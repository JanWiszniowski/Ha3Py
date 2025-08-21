"""
Module Bayesian by shift maximum magnitude computation
------------------------------------------------------

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
from ha3py.BayesianMmax import BayesianBase
from ha3py.ln_likelihood import ln_likelihood


class BayesianByShift(BayesianBase):
    r"""
    The class for estimation of the :math:`m_{max}` by the Bayesian by shift method.
    In the case of single catalog, the class define the likelihood function

    .. math::
        \mathcal{L}\left(\mathbf{m}|m_{max}\right)=
        \prod_{i=1}^{n}{f_M\left(m_i | m_{max}\right)},

    used in the Bayesian posterior distribution of :math:`m_{max}`

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

    and :math:`m_{max}^L={\widehat{m}}_{max}`.
    In the case of many catalogs :math:`\mathcal{L}` is computed by the ln_likelihood function.

    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict

    """
    def __init__(self, configuration, magnitude_distribution=None):
        """

        :param configuration:
        :type configuration:
        """
        self.configuration = configuration
        super().__init__(configuration, magnitude_distribution=magnitude_distribution)

    def _likelihood(self, m_max):
        """
        ln_likelihood(None, configuration, m_max=None)
        :param m_max:
        :type m_max:
        :return:
        :rtype:
        """
        if np.isscalar(m_max):
            return np.exp(-ln_likelihood(None, self.configuration, m_max=m_max))
        else:
            n = m_max.size
            output = []
            for idx in range(n):
                output.append(np.exp(ln_likelihood(None, self.configuration, m_max=m_max[idx])))
            return np.array(output)


def get_bayesian_by_shift(configuration, magnitude_distribution=None):
    """

    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict
    :param magnitude_distribution:
    :type magnitude_distribution:
    :return: The BayesianByShift object
    :rtype: BayesianBase derived object
    """
    return BayesianByShift(configuration, magnitude_distribution=magnitude_distribution)
