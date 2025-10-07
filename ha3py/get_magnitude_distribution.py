"""
Ha3Py
(c) Jan Wiszniowski, Andrzej Kijko
ver. 2024-01
"""

from ha3py.GutenbergRichter import GutenbergRichter
from ha3py.CompoundGutenbergRichter import CompoundGutenbergRichter
from ha3py.NonParametricPseudoGaussian import NonparametricGaussianKernel
from ha3py.utils import HaPyException


def get_magnitude_distribution(configuration, m_min=None, m_max=None, theta=None):
    r"""
    Function get_magnitude_distribution returns magnitude distribution object
    based on the parameter 'magnitude distribution' in the dictionary of Ha3Py Parameters
    There are three currently available classes of magnitude distribution objects:
    'Gutenberg-Richter'

    .. math::
        f_M(m)=\frac{\beta \exp[-\beta(m-m_{min})]}{1-\exp[-\beta(m_{max}-m_{min}])}

    .. math::
        F_M(m)=\frac{1-\exp[-\beta(m-m_{min}])}{1-\exp[-\beta(m_{max}-m_{min}])}

    'Compound Gutenberg-Richter'

    .. math::
        f_M(m) = \overline{\beta} C_{\beta}\left ( \frac{p}{p+m-m_{min}} \right )^{q+1}

    .. math::
        F_M(m)=C_{\beta} \left [1-\left ( \frac{p}{p+m-m_{min}} \right )^q  \right ]

    'Nonparametric gaussian kernel'

    .. math::
        f_M\left( m \right)=
        \frac{\left( h\sqrt{2\pi}\right)^{-1}\sum_{i=1}^{M}\frac{1}{T_i}\exp
        \left[ -0.5\left( \frac{m-m_i}{h} \right)^2 \right]}
        {\sum_{i=1}^{M}\frac{1}{T_i}\left[  \Phi\left( \frac{m_{max}-m_i}{h} \right)-
        \Phi\left( \frac{m_{min}-m_i}{h} \right) \right]}

    .. math::
        F_M\left( m \right)=
        \frac{\sum_{i=1}^{M}\frac{1}{T_i}\left[  \Phi\left( \frac{m-m_i}{h} \right)-
        \Phi\left( \frac{m_{min}-m_i}{h} \right) \right]}
        {\sum_{i=1}^{M}\frac{1}{T_i}\left[  \Phi\left( \frac{m_{max}-m_i}{h} \right)-
        \Phi\left( \frac{m_{min}-m_i}{h} \right) \right]}

    A more detailed description of the magnitude distribution is included in the definition of a specific class

    :param theta:
    :type theta:
    :param m_max:
    :param m_min:
    :param configuration: The dictionary of all Ha3Py parameters
    :type configuration: dict
    :return: the magnitude distribution object
    """
    beta = None
    if theta is not None and len(theta) >= 1:
        beta = theta[-1]
    distribution_name = configuration.get('magnitude_distribution', 'Nonparametric gaussian kernel')
    if distribution_name == 'Gutenberg-Richter':
        return GutenbergRichter(configuration, beta=beta, m_min=m_min, m_max=m_max)
    elif distribution_name == 'Compound Gutenberg-Richter':
        return CompoundGutenbergRichter(configuration, beta=beta, m_min=m_min, m_max=m_max)
    elif distribution_name == 'Nonparametric gaussian kernel':
        return NonparametricGaussianKernel(configuration, m_min=m_min, m_max=m_max)
    else:
        HaPyException('Unknown magnitude distribution')
