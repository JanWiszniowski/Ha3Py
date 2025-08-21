r"""
Return period calculation
-------------------------

:copyright:
    Jan Wiszniowski <jwisz@igf.edu.pl>,
    Andrzej Kijko <andrzej.kijko@up.ac.za>
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
:version 0.0.1:
    2025-01-01
"""


def return_period(m, lamb, mag_dist):
    r"""
    Procedure calculates the return period as the inverse
    of the annual probability of not occurrence,
    which is the survive function magnitude not occurrence distribution.

    .. math::
        \color{gray}T_{R}\left(m,\lambda_{0},F_{M}\right)=\frac{1}{S_M^{max}\left(m\right)}

    **Example**

    Let

    .. math::
        \color{gray}\lambda\left(m\right)=\lambda_{0}\left(1-F_{M}\left(m|m_{0}\right)\right),

    then return period is

    .. math::
        \color{gray}T_{R}\left(m,\lambda_{0},F_{M}\right)=\frac{1}{\lambda\left(m\right)}

    :param m: magnitude for the return period
    :param lamb: lambda for the minimum magnitude
    :param mag_dist: the magnitude distribution
    :return: the return period
    """
    lambda_x = lamb * mag_dist.sf(m)
    if lambda_x <= 1e-6:
        lambda_x = 1e-6
    return lambda_x, 1.0 / lambda_x


def grad_return_period(m, lamb, mag_distr):
    r"""
    Procedure calculates the returns gradient of return period parameters

    **Example**

    Gradient of return period for
    :math:`\color{gray}\lambda\left(m\right)=\lambda_{0}\left(1-F_{M}\left(m|m_{0}\right)\right)`:

    .. math::
        \color{gray}
        \frac{\partial T_{R}\left(m\right)}{\partial\lambda_{0}}=
        \frac{\partial\left(\frac{1}{\lambda_{0}\cdot\left(1-F_{M}\left(m|m_{0}\right)\right)}\right)}{\partial\lambda_{0}}=
        \frac{-1}{\lambda_{0}^{2}\cdot\left(1-F_{M}\left(m|m_{0}\right)\right)}

    Gradients of return period for typical parameters of a magnitude distribution.

    :math:`\beta`:

    .. math::
        \color{gray}
        \frac{\partial T_{R}\left(m\right)}{\partial\beta}=
        \frac{\partial\left(\frac{1}{\lambda_{0}\cdot\left(1-F_{M}\left(m|m_{0}\right)\right)}\right)}{\partial\beta}=
        \frac{-1}{\lambda\cdot\left(1-F_{M}\left(m|m_{0}\right)\right)^{2}}
        \frac{\partial\left(1-F_{M}\left(m|m_{0}\right)\right)}{\partial\beta}

    :math:`m_{max}`:

    .. math::
        \color{gray}
        \frac{\partial T_{R}\left(m\right)}{\partial m_{max}}=
        \frac{\partial\left(\frac{1}{\lambda_{0}\cdot\left(1-F_{M}\left(m|m_{0}\right)\right)}\right)}{\partial m_{max}}=
        \frac{-1}{\lambda\cdot\left(1-F_{M}\left(m|m_{0}\right)\right)^{2}}
        \frac{\partial\left(1-F_{M}\left(m|m_{0}\right)\right)}{\partial m_{max}}

    :param m: magnitude for the return period
    :param lamb: lambda for the minimum magnitude
    :param mag_distr: the magnitude distribution
    :return: gradient of the return priod versus lambad and magnitude distribution paramters
    """
    sf = mag_distr.sf(m)  # sf(m) = [1-cdf(m)]
    if sf == 0.0:
        print(f"Warning SF({m}) == 0")
        return None
    g_sf = mag_distr.grad_sf(m)  # grad sf = grad P(x>m)
    ret_dict = {'lambda': -1.0 / lamb / lamb / sf}
    # print(f"m{m}, lamb{lamb}, sf{sf}")
    for parx, gradx in g_sf.items():
        ret_dict[parx] = -1.0 / sf / sf / lamb * gradx
    return ret_dict
