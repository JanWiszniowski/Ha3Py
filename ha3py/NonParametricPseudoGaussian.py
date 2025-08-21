"""
Non-parametric magnitude distribution
-------------------------------------

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
from scipy import ndimage
from scipy.stats import norm
from ha3py.MagnitudeDistribution import BaseMagnitudeDistribution


def smooth_parameter_estimation(earthquakes):
    r"""
    Calculation of optimal smoothing parameter :math:`h` in the gaussian kernel function.

    Criterion of :math:`h` selection is `Adaptive Estimate of Spread`.

    .. math::
        h = 0.9 * \min\left[ sd,\frac{ Q\left( p=0.75\right)-Q\left( p=0.25\right)}{1.34}\right] n^{-1/5},

    where :math:`sd` is the standard deviation of magnitude values,
    :math:`Q\left( p=k\right)` is :math:`k` quantile,
    and :math:`n` is the number of considered events

    For more details, see Silverman (1986), pp. 48.

    :param earthquakes: column vector of earthquakes dictionaries
    :type earthquakes: list(dict)
    :return: optimal value of smoothing parameter in the non-parametric (Gaussian kernel) function
    :rtype: (float)
    """
    no_earthquakes = len(earthquakes)
    mag = [m['magnitude'] for m in earthquakes]
    mag.sort()
    i1 = round(0.25 * no_earthquakes)
    i2 = round(0.75 * no_earthquakes)
    quar_range = (mag[i2] - mag[i1]) / 1.34
    std_mag = ndimage.standard_deviation(mag)
    a = min(quar_range, std_mag)
    h = 0.9 * a * no_earthquakes ** (-0.2)
    # h = round(h, 2)
    return h


def msort(e):
    return e['magnitude']


class NonparametricGaussianKernel(BaseMagnitudeDistribution):
    r"""
    The continuous random variable.

    The non-parametric magnitude distribution with Gaussian kernel probability density function is:

    .. math::
        f_{M}(m)=\left\{ \begin{alignat*}{2}
        & 0 && : \text{for } m<m_{min} \\
        & \frac{\left( h\sqrt{2\pi}\right)^{-1}\sum_{i=1}^{M}\frac{1}{T_i}\exp
        \left[ -0.5\left( \frac{m-m_i}{h} \right)^2 \right]}
        {\sum_{i=1}^{M}\frac{1}{T_i}\left[  \Phi\left( \frac{m_{max}-m_i}{h} \right)-
        \Phi\left( \frac{m_{min}-m_i}{h} \right) \right]} &&
        : \text{for } m_{min} \leqslant m \leqslant m_{max}\\
        & 0 && : \text{for } m>m_{max}
        \end{alignat*} \right.

    where :math:`h` is the smoothing parameter and :math:`T_i` is the sum of time periods
    of catalogues, where :math:`m_i \ge m_x`. Values :math:`1/T_i` can be replaced by predefined weights
    if they are defined in catalogs. The smoothing parameter :math:`h` is estimated by
    the smooth_parameter_estimation function unless it is predefined.


    The non-parametric magnitude distribution with Gaussian kernel cumulate density functions:

    .. math::
        F_{M}(m)=\left\{ \begin{alignat*}{2}
        & 0 && : \text{for } m<m_{min} \\
        & \frac{\sum_{i=1}^{M}\frac{1}{T_i}\left[  \Phi\left( \frac{m-m_i}{h} \right)-
        \Phi\left( \frac{m_{min}-m_i}{h} \right) \right]}
        {\sum_{i=1}^{M}\frac{1}{T_i}\left[  \Phi\left( \frac{m_{max}-m_i}{h} \right)-
        \Phi\left( \frac{m_{min}-m_i}{h} \right) \right]} &&
        : \text{for } m_{min} \leqslant m \leqslant m_{max}\\
        & 1 && : \text{for } m>m_{max}
        \end{alignat*} \right.

    where :math:`\Phi(m)` is the normal cdf.

    The only variable parameter is :math:`m_{max}`. The survive function gradient for :math:`m_{max}` is

    .. math::
        \frac{\partial S_{M}(m)}{\partial m_{max}}=
        \frac{F_{M}(m)f_M\left( m_{max} \right)}{\sum_{i=1}^{M}\frac{1}{T_i}
        \left[  \Phi\left( \frac{m_{max}-m_i}{h} \right)-
        \Phi\left( \frac{m_{min}-m_i}{h} \right) \right]}

    """

    # def __init__(self, magnitudes, smooth_par=np.NaN, m_min=np.NaN, M_max=10):
    def __init__(self, configuration, no_largest=None, m_min=None, m_max=None):
        r"""
        Required Parameters:
            paleo_catalog
            historic_catalog
            complete_catalogs
        Optional Parameters:
            smoothing_factor
        """
        p_phs = configuration.get('paleo_catalog')
        p_his = configuration.get('historic_catalog')
        p_comp = configuration.get('complete_catalogs')
        no_largest = configuration.get('no_largest_magnitudes', no_largest)
        mag_v = []
        catalog_parameters = []
        if p_phs:
            for m in p_phs['earthquakes']:
                mag_v.append(m)
            catalog_parameters.append((p_phs['m_min'], p_phs['time_span']))
        if p_his:
            for m in p_his['earthquakes']:
                mag_v.append(m)
            catalog_parameters.append((p_his['m_min'], p_his['time_span']))
        if p_comp:
            for cat in p_comp:
                for m in cat['earthquakes']:
                    mag_v.append(m)
                catalog_parameters.append((cat['m_min'], cat['time_span']))
        mag_v.sort(key=msort, reverse=True)
        if no_largest is not None and no_largest < len(mag_v):
            mag_v = mag_v[0:no_largest]
        mag = list()
        for v in mag_v:
            m = v['magnitude']
            w = v.get('weight')
            if w is None:
                t = 0
                for event_parameters in catalog_parameters:
                    if m >= event_parameters[0]:
                        t += event_parameters[1]
                w = 1.0 / t
            mag.append((m, w))
        self.mag = mag
        self.n = len(mag)
        if self.n < 2:
            raise Exception("To few magnitudes")
        self.h = configuration.get('smoothing_factor', smooth_parameter_estimation(mag_v))
        # self.p_min = None
        # self.den = None
        super().__init__(configuration, 'N-P-G-K', long_name='Gaussian kernel non-parametric magnitude distribution',
                         m_min=m_min, m_max=m_max)

    def _prepare(self):
        p_min = 0.0
        p_max = 0.0
        sig = self.h
        for x_i in self.mag:
            # p_min += f_cdf_gauss2(self.m_min, x_i[0], sig) * x_i[1]
            # p_max += f_cdf_gauss2(self.m_max, x_i[0], sig) * x_i[1]
            p_min += norm.cdf(self._m_min, x_i[0], sig) * x_i[1]
            p_max += norm.cdf(self._m_max, x_i[0], sig) * x_i[1]
        self.p_min = p_min
        self.den = p_max - p_min

    def _cdf(self, x, *args):
        p = np.zeros(len(x), dtype=float)
        sig = self.h
        # sig5 = sig * 5.0
        for x_i in self.mag:
            # dx = x - x_i[0]
            # if dx > sig5:
            #     p += 1.0
            # elif dx > -sig5:
            #     p = p + norm.cdf(x, x_i[0], sig) * x_i[1]
            p = p + norm.cdf(x, x_i[0], sig) * x_i[1]
        return (p-self.p_min) / self.den

    def _pdf(self, x, *args):
        p = np.zeros(len(x), dtype=float)
        sig = self.h
        # sig5 = sig * 5.0
        for x_i in self.mag:
            # if np.abs(x - x_i[0]) < sig5:
            #     p = p + norm.pdf(x, x_i[0], sig) * x_i[1]
            p = p + norm.pdf(x, x_i[0], sig) * x_i[1]
        return p / self.den

    def _grad_sf(self, m):
        grad = self._cdf(m) / self.den * self.pdf(self.m_max)
        return {'m_max': grad}

    def _coefficient_names(self):
        return ['m_max']

    def _coefficient_values(self):
        return [self.m_max]

    def _const_coefficients(self):
        return ['h']

class NonparametricPseudoGaussianKernel(BaseMagnitudeDistribution):
    r"""
    The continuous random variable.

    The non-parametric magnitude distribution with pseudo Gaussian kernel probability density function is:

    .. math::
        f_{M}(m)=\left\{ \begin{alignat*}{2}
        & 0 && : \text{for } m<m_{min} \\
        & \frac{\sum_{i=1}^{M}\frac{1}{T_i} \psi\left(\frac{m-m_i}{h}\right)}
        {\sum_{i=1}^{M}\frac{1}{T_i}\left[  \Psi\left( \frac{m_{max}-m_i}{h} \right)-
        \Psi\left( \frac{m_{min}-m_i}{h} \right) \right]} &&
        : \text{for } m_{min} \leqslant m \leqslant m_{max}\\
        & 0 && : \text{for } m>m_{max}
        \end{alignat*} \right.

    where

    .. math::
        \psi(x) = 2\frac{(a_1 + 2a_2x + 3a_3x^2 + 4a_4x^3)}{ \sigma(1 + a_1x + a_2x^2 + a_3x^3 + a_4x^4)^5},

        \Psi(x) = 1.0 - 0.5(1 + a_1x + a_2x^2 + a_3x^3 + a_4x^4)^{-4}

        a_1 = 0.196854,
        a_2 = 0.115194,
        a_3 = 0.000344,
        a_4 = 0.019527,

    :math:`h` is the smoothing parameter and :math:`T_i` is the sum of time periods
    of catalogues, where :math:`m_i \ge m_x`. Values :math:`1/T_i` can be replaced by predefined weights
    if they are defined in catalogs. The smoothing parameter :math:`h` is estimated by
    the smooth_parameter_estimation function unless it is predefined.


    The non-parametric magnitude distribution with Gaussian kernel cumulate density functions:

    .. math::
        F_{M}(m)=\left\{ \begin{alignat*}{2}
        & 0 && : \text{for } m<m_{min} \\
        & \frac{\sum_{i=1}^{M}\frac{1}{T_i}\left[  \Psi\left( \frac{m-m_i}{h} \right)-
        \Psi\left( \frac{m_{min}-m_i}{h} \right) \right]}
        {\sum_{i=1}^{M}\frac{1}{T_i}\left[  \Psi\left( \frac{m_{max}-m_i}{h} \right)-
        \Psi\left( \frac{m_{min}-m_i}{h} \right) \right]} &&
        : \text{for } m_{min} \leqslant m \leqslant m_{max}\\
        & 1 && : \text{for } m>m_{max}
        \end{alignat*} \right.

    The only variable parameter is :math:`m_{max}`. The survive function gradient for :math:`m_{max}` is

    .. math::
        \frac{\partial S_{M}(m)}{\partial m_{max}}=
        \frac{F_{M}(m)f_M\left( m_{max} \right)}{\sum_{i=1}^{M}\frac{1}{T_i}
        \left[  \Psi\left( \frac{m_{max}-m_i}{h} \right)-
        \Psi\left( \frac{m_{min}-m_i}{h} \right) \right]}

    """

    # def __init__(self, magnitudes, smooth_par=np.NaN, m_min=np.NaN, M_max=10):
    def __init__(self, configuration, no_largest=None, m_min=None, m_max=None):
        r"""
        Required Parameters:
            paleo_catalog
            historic_catalog
            complete_catalogs
        Optional Parameters:
            smoothing_factor
        """
        p_phs = configuration.get('paleo_catalog')
        p_his = configuration.get('historic_catalog')
        p_comp = configuration.get('complete_catalogs')
        no_largest = configuration.get('no_largest_magnitudes', no_largest)
        mag_v = []
        catalog_parameters = []
        if p_phs:
            for m in p_phs['earthquakes']:
                mag_v.append(m)
            catalog_parameters.append((p_phs['m_min'], p_phs['time_span']))
        if p_his:
            for m in p_his['earthquakes']:
                mag_v.append(m)
            catalog_parameters.append((p_his['m_min'], p_his['time_span']))
        if p_comp:
            for cat in p_comp:
                for m in cat['earthquakes']:
                    mag_v.append(m)
                catalog_parameters.append((cat['m_min'], cat['time_span']))
        mag_v.sort(key=msort, reverse=True)
        if no_largest is not None and no_largest < len(mag_v):
            mag_v = mag_v[0:no_largest]
        mag = list()
        for v in mag_v:
            m = v['magnitude']
            w = v.get('weight')
            if w is None:
                t = 0
                for event_parameters in catalog_parameters:
                    if m >= event_parameters[0]:
                        t += event_parameters[1]
                w = 1.0 / t
            mag.append((m, w))
        self.mag = mag
        self.n = len(mag)
        if self.n < 2:
            raise Exception("To few magnitudes")
        self.h = configuration.get('smoothing_factor', smooth_parameter_estimation(mag_v))
        # self.p_min = None
        # self.den = None
        super().__init__(configuration, 'N-P-G-K', long_name='Gaussian kernel non-parametric magnitude distribution',
                         m_min=m_min, m_max=m_max)

    def _prepare(self):
        p_min = 0.0
        p_max = 0.0
        sig = self.h
        for x_i in self.mag:
            # p_min += f_cdf_gauss2(self.m_min, x_i[0], sig) * x_i[1]
            # p_max += f_cdf_gauss2(self.m_max, x_i[0], sig) * x_i[1]
            p_min += norm.cdf(self._m_min, x_i[0], sig) * x_i[1]
            p_max += norm.cdf(self._m_max, x_i[0], sig) * x_i[1]
        self.p_min = p_min
        self.den = p_max - p_min

    def _cdf(self, x, *args):
        a1 = .196854
        a2 = .115194
        a3 = .000344
        a4 = .019527

        p = np.zeros(len(x), dtype=float)
        sig = self.h
        for x_i in self.mag:
            xn = (x - x_i[0]) / sig
            xa = np.abs(xn)

            xa2 = np.power(xa,2)
            xa3 = np.power(xa,3)
            xa4 = np.power(xa,4)

            y = 1.0 - 0.5 * (1 + a1 * xa + a2 * xa2 + a3 * xa3 + a4 * xa4)**(-4)
            p = p + np.where(xn < 0.0, 1.0-y, y) * x_i[1]
        return (p-self.p_min) / self.den

    def _pdf(self, x, *args):
        a1 = .196854
        a2 = .115194
        a3 = .000344
        a4 = .019527

        p = np.zeros(len(x), dtype=float)
        sig = self.h
        for x_i in self.mag:
            xn = (x - x_i[0]) / sig
            xa = np.abs(xn)

            xa2 = np.power(xa,2)
            xa3 = np.power(xa,3)
            xa4 = np.power(xa,4)

            y = 2.0 * (a1 + 2.0*a2*xa + 3.0*a3*xa2 + 4.0*a4*xa3)*(1 + a1*xa + a2*xa2 + a3*xa3 + a4*xa4)**(-5) / sig
            p = p + np.where(xn < 0.0, -y, y) * x_i[1]
        return p / self.den

    def _grad_sf(self, m):
        grad = self._cdf(m) / self.den * self.pdf(self.m_max)
        return {'m_max': grad}

    def _coefficient_names(self):
        return ['m_max']

    def _coefficient_values(self):
        return [self.m_max]

    def _const_coefficients(self):
        return ['h']