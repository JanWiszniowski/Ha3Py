"""
Module Compound Gutenberg-Richter magnitude distribution
--------------------------------------------------------

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

from ha3py.MagnitudeDistribution import BaseMagnitudeDistribution


class CompoundGutenbergRichter(BaseMagnitudeDistribution):
    r"""
    The continuous random variable.

    The compound Gutenberg-Richter-Bayes probability density function is:

    .. math::
        f_{M}(m)=\left\{ \begin{alignat*}{2}
        & 0 && : \text{for } m<m_{min} \\
        & \overline{\beta} C_{\beta}\left[ \frac{q_{\beta}}{q_{\beta}+
        \overline{\beta}\left(m-m_{min} \right)}\right ]^{q_{\beta}+1} &&
        : \text{for } m_{min} \leqslant m \leqslant m_{max}\\
        & 0 && : \text{for } m>m_{max}
        \end{alignat*} \right.

    where:

    .. math::
        C_{\beta}=\frac{1}
        {1- \left [\frac{q_{\beta}}{q_{\beta}+\overline{\beta}\left(m_{max}-m_{min}\right)} \right ]^{q_{\beta}}}

    and :math:`q_{\beta}=\left ( \overline{\beta}/\sigma _{\beta} \right )^2=\text{constant}`.

    The compound Gutenberg-Richter cumulate density function is:

    .. math::
        F_{M}(m)=\left\{ \begin{alignat*}{2}
        & 0 && : \text{for } m<m_{min} \\
        & C_{\beta} \left [1-\left ( \frac{q_{\beta}}{q_{\beta}+\overline{\beta}\left(m-m_{min}\right)} \right )
        ^{q_{\beta}}  \right ] &&
        : \text{for } m_{min} \leqslant m \leqslant m_{max}\\
        & 1 && : \text{for } m>m_{max}
        \end{alignat*} \right.

    Gradient of the survive function for :math:`m_{max}` is:

    .. math::
        \frac{\partial(1-F_{M}\left(m\right))}{\partial m_{max}}=
        -\frac{N}{D^{2}}\overline{\beta}\left(\frac{q_{\beta}}
        {q_{\beta}+\overline{\beta}\left(m_{max}-m_{min}\right)}\right)^{q_{\beta}+1}

    Gradient of the survive function for :math:`\beta` is

    .. math::
        \frac{\partial(1-F_{M}(m))}{\partial\overline{\beta}}=\frac
        {N\left(m_{max}-m_{min}\right)\left(\frac{q_{\beta}}{q_{\beta}+\overline{\beta}\left(m_{max}
        -m_{min}\right)}\right)^{q_{\beta}+1}-D\left(m-m_{min}\right)\left(\frac{q_{\beta}}
        {q_{\beta}+\overline{\beta}\left(m-m_{min}\right)}\right)^{q_{\beta}+1}}{D^{2}}

    """
    # def __init__(self, _beta=np.NaN, b_std=np.NaN, p=np.NaN, q=np.NaN, m_min=0, M_max=10):
    def __init__(self, parameters, beta=None, q_beta=None, m_min=None, m_max=None):
        r"""
        Required params:
            m_min
            m_max_current
            _beta
        Optional gr_parameters:
            M_max (if m_max_current missing)
            sd_beta (if q_beta missing)
        """
        if beta is not None:
            self.beta = beta
        else:
            self.beta = parameters.get('beta', 2.3)
        if q_beta is not None:
            self.q = q_beta
        else:
            self.q = parameters.get('q_beta', 100.0)
        super().__init__(parameters, 'C-G-R', 'Compound Gutenberg-Richter magnitude distribution',
                         m_min=m_min, m_max=m_max)

    def _prepare(self):
        self.den = (1 - (self.q/(self.q + self.beta * (self.m_max - self.m_min)))**self.q)

    def _pdf(self, m, *args):
        r"""

        .. math:
            \text{pdf}\left( m \right) =
            \beta\left[ \frac{q_\beta}{q_\beta+\beta\left( m-m_{min}\right) } \right]^{q_\beta+1}/D_\beta

        where

        .. math:
            D_\beta =
            \left[ 1-\frac{q_\beta}{q_\beta+\beta\left( m_{max}-m_{min}\right) } \right]^{q_\beta}

        :param m:
        :param args:
        :return:
        """
        return self.beta * (self.q / (self.q + self.beta * (m - self.m_min)))**(self.q + 1)/self.den

    def _cdf(self, m, *args):
        r"""

        .. math:
        \text{cdf}\left( m \right) = \left[ 1-\frac{q_\beta}{q_\beta+\beta\left( m-m_{min}\right) } \right]^
        {q_\beta}/D_\beta
        where
        .. math:
        D_\beta = \left[ 1-\frac{q_\beta}{q_\beta+\beta\left( m_{max}-m_{min}\right) } \right]^{q_\beta}

        :param m:
        :param args:
        :return:
        """
        return (1 - (self.q / (self.q + self.beta * (m - self.m_min))) ** self.q) / self.den

    def _grad_sf(self, m):
        r"""
        Gradient of the survive function for :math:`m_{max}` is:

        .. math::
            \frac{\partial(1-F_{M}\left(m\right))}{\partial m_{max}}=
            -\frac{N}{D^{2}}\overline{\beta}\left(\frac{q_{\beta}}
            {q_{\beta}+\overline{\beta}\left(m_{max}-m_{min}\right)}\right)^{q_{\beta}+1}

        Gradient of the survive function for :math:`\beta` is

        .. math::
            \frac{\partial(1-F_{M}(m))}{\partial\overline{\beta}}=\frac
            {N\left(m_{max}-m_{min}\right)\left(\frac{q_{\beta}}{q_{\beta}+\overline{\beta}\left(m_{max}
            -m_{min}\right)}\right)^{q_{\beta}+1}-D\left(m-m_{min}\right)\left(\frac{q_{\beta}}
            {q_{\beta}+\overline{\beta}\left(m-m_{min}\right)}\right)^{q_{\beta}+1}}{D^{2}}

        :param m:
        :return:
        """
        if self.m_max >= m >= self.m_min:
            # nominator
            diff_m = m - self.m_min
            nom_pow1 = (self.q/(self.q + self.beta*diff_m))**(self.q+1)
            nom = 1 - (self.q/(self.q + self.beta*diff_m))**self.q
            # denominatar
            diff_m_max = self.m_max - self.m_min
            den = self.den
            den_pow1 = (self.q/(self.q + self.beta*diff_m_max))**(self.q+1)
            # _beta
            d_beta = (nom * diff_m_max * den_pow1 - den * diff_m * nom_pow1) / den / den
            # M_max
            d_m_max = nom * self.beta * den_pow1 / den / den
        else:
            d_beta = 0
            d_m_max = 0
        return {'beta': d_beta, 'm_max': d_m_max}

    def _coefficient_names(self):
        return ['beta', 'm_max']

    def _coefficient_values(self):
        return [self.beta, self.m_max]

    def _const_coefficients(self):
        return ['q_beta', 'm_min']


if __name__ == "__main__":
    Parameters = {'beta': 2.303, 'q_beta': 10, 'm_min': 2.0, 'M_max': 6.0}
    mag_distr = CompoundGutenbergRichter(Parameters)
    magnitudes = [x / 10 for x in range(0, 71)]
    # print('magnitudes = {}'.format(m))
    # cdf_vals =  magnitude_distribution.cdf(m)
    # print('cdf = {}'.format(cdf_vals))
    # pdf_vals =  magnitude_distribution.pdf(m)
    # print('pdf = {}'.format(pdf_vals))
    for mag in magnitudes:
        print('magnitudes = {}'.format(mag))
        cdf_vals = mag_distr.cdf(mag)
        print('cdf = {}'.format(cdf_vals))
        pdf_vals = mag_distr.pdf(mag)
        print('pdf = {}'.format(pdf_vals))
    # mag = 4
    # print('magnitudes = {}'.format(mag))
    # cdf_vals = magnitude_distribution.cdf(mag)
    # print('cdf = {}'.format(cdf_vals))
    # pdf_vals = magnitude_distribution.pdf(mag)
    # print('pdf = {}'.format(pdf_vals))
