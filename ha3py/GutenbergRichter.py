"""
Gutenberg-Richter magnitude distribution
----------------------------------------

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
from ha3py.MagnitudeDistribution import BaseMagnitudeDistribution


class GutenbergRichter(BaseMagnitudeDistribution):
    r"""
    The continuous random variable.

    The Gutenberg-Richter probability density function is:

    .. math::
        f_{M}(m)=\begin{cases}
        & 0 && : \text{for } m<m_{min} \\
        & \frac {\beta\exp\left[-\beta\left(m-m_{min}\right)\right]}
        {1-\exp\left[-\beta\left(m_{max}-m_{min}\right)\right]} &&
        : \text{for } m_{min} \leqslant m \leqslant m_{max}\\
        & 0 && : \text{for } m>m_{max}
        \end{cases}

    The Gutenberg-Richter cumulate density function is

    .. math::
        F_{M}(m)=\left\{ \begin{alignat*}{2}
        & 0 && : \text{for } m<m_{min} \\
        & \frac{1-\exp\left[-\beta\left(m-m_{min}\right)\right]}
        {1-\exp\left[-\beta\left(m_{max}-m_{min}\right)\right]} &&
        : \text{for } m_{min} \leqslant m \leqslant m_{max}\\
        & 1 && : \text{for } m>m_{max}
        \end{alignat*} \right.

    The Gutenberg-Richter gradients for :math:`m_{max}` and :math:`\beta` are:

    .. math::
        \frac{\partial S_{M}(m)}{\partial m_{max}}=\frac
        {N\beta\exp\left[-\beta\left(m_{max}-m_{min}\right)\right]} {D^{2}}

    and

    .. math::
        \frac{\partial S_{M}(m)}{\partial\beta}=\frac
        {N\left(m_{max}-m_{min}\right)\exp\left[-\beta\left(m_{max}-m_{min}\right)\right]
        -D\left(m-m_{min}\right)\exp\left[-\beta\left(m-m_{min}\right)\right]} {D^{2}}

    where

    .. math::
        D=1-\exp\left[-\beta\left(m_{max}-m_{min}\right)\right]

    end

    .. math::
        N=1-\exp\left[-\beta\left(m-m_{min}\right)\right]

    """

    def __init__(self, parameters, beta=None, m_min=None, m_max=None):
        self.den = 1.0
        if beta:
            self._beta = beta
        else:
            self._beta = parameters.get('beta', 2.3)
        super().__init__(parameters, 'G-R', long_name='Gutenberg-Richter magnitude distribution',
                         m_min=m_min, m_max=m_max)

    def _prepare(self):
        self.den = (1 - np.exp(-self._beta * (self.m_max - self.m_min)))

    def _pdf(self, x, *args):
        r"""
        Returns the probability density function

        .. math::
            f_{M}(m)=\left\{ \begin{alignat*}{2}
            & 0 && : \text{for } m<m_{min} \\
            & \frac {\beta\exp\left[-\beta\left(m-m_{min}\right)\right]}
            {1-\exp\left[-\beta\left(m_{max}-m_{min}\right)\right]} &&
            : \text{for } m_{min} \leqslant m \leqslant m_{max}\\
            & 0 && : \text{for } m>m_{max}
            \end{alignat*} \right.

        """
        return self._beta * np.exp(-self._beta * (x - self.m_min)) / self.den

    def _cdf(self, x, *args):
        r"""
        Returns the cumulative density function
        .. math::
            F_{M}(m)=\left\{ \begin{alignat*}{2}
            & 0 && : \text{for } m<m_{min} \\
            & \frac{1-\exp\left[-\beta\left(m-m_{min}\right)\right]}
            {1-\exp\left[-\beta\left(m_{max}-m_{min}\right)\right]} &&
            : \text{for } m_{min} \leqslant m \leqslant m_{max}\\
            & 1 && : \text{for } m>m_{max}
            \end{alignat*} \right.

        """
        return (1 - np.exp(-self._beta * (x - self.m_min))) / self.den

    def _grad_sf(self, m):
        r"""
        The Gutenberg-Richter gradients for :math:`m_{max}` and :math:`\beta` are:

        .. math::
            \frac{\partial S_{M}(m)}{\partial m_{max}}=\frac
            {N\beta\exp\left[-\beta\left(m_{max}-m_{min}\right)\right]} {D^{2}}

        and

        .. math::
            \frac{\partial S_{M}(m)}{\partial\beta}=\frac
            {N\left(m_{max}-m_{min}\right)\exp\left[-\beta\left(m_{max}-m_{min}\right)\right]
            -D\left(m-m_{min}\right)\exp\left[-\beta\left(m-m_{min}\right)\right]} {D^{2}}

        where

        .. math::
            D=1-\exp\left[-\beta\left(m_{max}-m_{min}\right)\right]

        end

        .. math::
            N=1-\exp\left[-\beta\left(m-m_{min}\right)\right]

        :param m:
        :return:
        """
        if self.m_max >= m >= self.m_min:
            # nominator
            diff_m = m - self.m_min
            nom_exp = np.exp(-self._beta * diff_m)
            nom = 1 - nom_exp
            # denominator
            diff_m_max = self.m_max - self.m_min
            den_exp = np.exp(-self._beta * diff_m_max)
            den = self.den
            # _beta
            d_beta = (nom * diff_m_max * den_exp - den * diff_m * nom_exp) / den / den
            # M_max
            d_m_max = nom * self._beta * den_exp / den / den
        else:
            d_beta = 0
            d_m_max = 0
        return {'beta': d_beta, 'm_max': d_m_max}

    def _coefficient_names(self):
        return ['beta', 'm_max']

    def _coefficient_values(self):
        return [self.beta, self.m_max]

    def _const_coefficients(self):
        return ['m_min']

    @property
    def beta(self):
        """It is the G-R beta value"""
        return self._beta

    @beta.setter
    def beta(self, val):
        self._beta = val
        self._prepare()

    @beta.getter
    def beta(self):
        return self._beta


if __name__ == "__main__":
    gr_parameters = {'beta': 2.303, 'q_beta': 10, 'm_min': 2.0, 'm_max': 6.0}
    magnitude_distribution = GutenbergRichter(gr_parameters)
    magnitudes = np.array([x / 10 for x in range(0, 71)])
    mag = 4
    print(f'magnitudes = {mag}')
    cdf_vals = magnitude_distribution.cdf(mag)
    print(f'cdf = {cdf_vals}')
    pdf_vals = magnitude_distribution.pdf(mag)
    print(f'pdf = {pdf_vals}')
