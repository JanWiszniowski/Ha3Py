r"""
Plotting results procedures
---------------------------

:copyright:
    Jan Wiszniowski <jwisz@igf.edu.pl>,
    Andrzej Kijko <andrzej.kijko@up.ac.za>
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
:version 0.0.1:
    2025-01-01

The module contain procedure for plotting seismic hazard diagrams resulted
from estimation of earthquake hazard parameters by Ha3Py.
Most functions plot diagrams curves with uncertainty margins.
"""

import numpy as np
import matplotlib.pyplot as plt
from ha3py.return_period import return_period, grad_return_period
# from ha3py.define_prob import get_max_prob
from ha3py.configuration import load_configuration
from ha3py.get_events_occurrence import get_events_occurrence


def fill_sd(x, y, ax):
    # """
    # Procedure fills the area of uncertainty
    # :param x:
    # :type x:
    # :param y:
    # :type y:
    # :param ax:
    # :type ax:
    # :return:
    # :rtype:
    # """
    sd_color = (0.4, 0.3, 0.5)
    sd_facecolor = sd_color
    sd_alpha = 0.5
    ax.fill(x, y, color=sd_color, facecolor=sd_facecolor, alpha=sd_alpha)


def plot_mean(x, y, ax):
    # """
    # Procedure plit the mean value curve
    # :param x:
    # :type x:
    # :param y:
    # :type y:
    # :param ax:
    # :type ax:
    # :return:
    # :rtype:
    # """
    m_color = 'r'
    m_width = 2
    ax.plot(x, y, color=m_color, linewidth=m_width)


def comp_mean_return_period(configuration):
    """
    Procedure compute the mean return period

    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict
    :return:
    :rtype:
    """
    event_occurrence = get_events_occurrence(configuration)
    lamb = configuration['lambda']
    m_min = configuration['m_min']
    m_max = configuration['m_max']
    rp_max = 1.0e5
    d_mag = 0.005
    rp_v = []
    temp_m_max = round(m_max, 3)
    event_occurrence.m_max = temp_m_max
    m_v = np.arange(m_min, temp_m_max, d_mag)
    for mag in np.nditer(m_v):
        [_, rp] = return_period(mag, lamb, event_occurrence.magnitude_distribution)
        if rp_max >= rp >= 0:
            rp_v.append((mag, rp))
    return rp_v


def comp_sd_return_period(configuration, rp_sign):
    event_occurrence = get_events_occurrence(configuration)
    lamb = configuration['lambda']
    m_min = configuration['m_min']
    m_max = configuration['m_max']
    sd_m_max = configuration['sd_m_max']
    cov_mc = np.array(configuration['cov_beta_lambda'])
    rp_max = 1.0e5
    d_mag = 0.005
    rp_v = []
    last_rp_sd = 0.0
    temp_m_max = round(m_max - rp_sign * sd_m_max, 3)
    event_occurrence.m_max = temp_m_max
    m_v = np.arange(m_min, temp_m_max, d_mag)
    # m_v = np.arange(m_min, temp_m_max+d_mag/2, d_mag)
    coefficient_names = event_occurrence.coefficient_names[:-1]
    for mag in np.nditer(m_v):  # MAGNITUDE LOOP STARTS HERE ----------------
        [_, rp] = return_period(mag, lamb, event_occurrence.magnitude_distribution)
        # ----------------- Calculation of SD of Return Period --------------------
        grad_dict = grad_return_period(mag, lamb, event_occurrence.magnitude_distribution)
        if not grad_dict:
            continue
        # grad_vec = np.array([grad_dict['beta'], grad_dict['lambda']])[np.newaxis]
        grad_vec = np.array([grad_dict[name] for name in coefficient_names])[np.newaxis]
        sd_rp = np.sqrt(abs(np.matmul(grad_vec, np.matmul(cov_mc, grad_vec.T))))
        # -------------- END: Calculation of SD of Return period ------------------
        rp_sd = rp + sd_rp * rp_sign
        if rp_max >= rp_sd >= 0 and rp_sd >= last_rp_sd:
            last_rp_sd = rp_sd
            rp_v.append((mag, rp_sd[0][0]))
            # print(f"{mag}, {rp_sd[0][0]}")
    return rp_v


def plot_return_period(configuration, ax=None):
    rp_v = comp_mean_return_period(configuration)
    rp_psd_v = comp_sd_return_period(configuration, +1)
    rp_msd_v = comp_sd_return_period(configuration, -1)
    rp_max_max = max(rp_v[-1][1], rp_msd_v[-1][1], rp_psd_v[-1][1])
    rp_v.append((rp_v[-1][0], rp_max_max))
    rp_msd_v.append((rp_psd_v[-1][0], rp_max_max))
    rp_psd_v.append((rp_msd_v[-1][0], rp_max_max))
    sd_x, sd_y = [], []
    for x, y in rp_msd_v:
        sd_x.append(x)
        sd_y.append(y)
    rp_psd_v.reverse()
    for x, y in rp_psd_v:
        sd_x.append(x)
        sd_y.append(y)
    m_x, m_y = [], []
    for x, y in rp_v:
        m_x.append(x)
        m_y.append(y)
    if not ax:
        fig, ax = plt.subplots()
    fill_sd(sd_x, sd_y, ax)
    plot_mean(m_x, m_y, ax)
    ax.set_xlabel('Magnitude')
    ax.set_ylabel('Return period [yrs]')
    ax.set_title('Area: {}'.format(configuration['area_name']))
    ax.grid(which='both', linestyle=':')

    return ax


def comp_mean_prob(configuration, tp):
    event_occurrence = get_events_occurrence(configuration)
    m_min = configuration['m_min']
    m_max = configuration['m_max']
    cpr_min = 1e-6
    d_mag = 0.005
    cpr_v = []
    temp_m_max = round(m_max, 3)
    event_occurrence.m_max = temp_m_max
    m_v = np.arange(m_min, temp_m_max + d_mag / 2, d_mag)
    for mag in np.nditer(m_v):
        cpr = event_occurrence.sf(mag, tp)  # 1 - PROB.
        if cpr >= cpr_min:
            cpr_v.append((mag, cpr))
        # if round(float(mag), 1) == round(float(mag), 4):
        #     print(f"{mag:.3f}, {cpr:.5f}")
    return cpr_v


def comp_sd_prob(pars, tp, cpr_sign):
    event_occurrence = get_events_occurrence(pars)
    m_min = pars['m_min']
    m_max = pars['m_max']
    sd_m_max = pars['sd_m_max']
    cov_mc = np.array(pars['cov_beta_lambda'])
    cpr_min = 1e-6
    d_mag = 0.005
    cpr_v = []
    temp_m_max = round(m_max + cpr_sign * sd_m_max, 3)
    event_occurrence.m_max = temp_m_max
    m_v = np.arange(m_min, temp_m_max + d_mag / 2, d_mag)
    coefficient_names = event_occurrence.coefficient_names[:-1]
    for mag in np.nditer(m_v):
        cpr = event_occurrence.sf(mag, tp)  # 1 - PROB.
        grad_dict = event_occurrence.grad_sf(mag, 1.0)
        # grad_vec = np.array([grad_dict['beta'], grad_dict['lambda']])[np.newaxis]
        grad_vec = np.array([grad_dict[name] for name in coefficient_names])[np.newaxis]
        sd_cpr = np.sqrt(abs(np.matmul(grad_vec, np.matmul(cov_mc, grad_vec.T))))[0][0]
        cpr_sd = cpr + sd_cpr * cpr_sign
        # cpr_sd = 0.5 + sd_cpr * cpr_sign
        if cpr_sd > 1.0:
            cpr_sd = 1.0
        if cpr_sd >= cpr_min:
            cpr_v.append((mag, cpr_sd))
        # if round(float(mag), 1) == round(float(mag), 4):
        #     print(f"{mag:.3f}, {cpr_sd:.5f}")
    return cpr_v


def plot_prob(pars, tp, ax=None):
    cpr_v = comp_mean_prob(pars, tp)
    cpr_psd_v = comp_sd_prob(pars, tp, +1)
    cpr_msd_v = comp_sd_prob(pars, tp, -1)
    # m1 = min(cpr_v, key=lambda tup: tup[1])
    # m2 = min(cpr_msd_v, key=lambda tup: tup[1])
    # m3 = min(cpr_psd_v, key=lambda tup: tup[1])
    cpr_min_min = min(cpr_v[-1][1], cpr_msd_v[-1][1], cpr_psd_v[-1][1])
    # cpr_min_min = min(min(cpr_v)[1], min(cpr_msd_v)[1], min(cpr_psd_v)[1])
    cpr_v.append((cpr_v[-1][0], cpr_min_min))
    cpr_msd_v.append((cpr_msd_v[-1][0], cpr_min_min))
    cpr_psd_v.append((cpr_psd_v[-1][0], cpr_min_min))
    sd_x, sd_y = [], []
    for x, y in cpr_msd_v:
        sd_x.append(x)
        sd_y.append(y)
    cpr_psd_v.reverse()
    for x, y in cpr_psd_v:
        sd_x.append(x)
        sd_y.append(y)
    m_x, m_y = [], []
    for x, y in cpr_v:
        m_x.append(x)
        m_y.append(y)
    if not ax:
        fig, ax = plt.subplots()
    fill_sd(sd_x, sd_y, ax)
    plot_mean(m_x, m_y, ax)
    ax.set_xlabel('Magnitude')
    ax.set_ylabel('Annual probability of exceedance')
    ax.set_title('Area: {}'.format(pars['area_name']))
    ax.grid(which='both', linestyle=':')
    return ax


def plot_hazard(pars, what, ax=None):
    # if 'hazard' not in pars:
    #     return
    if not ax:
        fig, ax = plt.subplots()
    # mags = [item['mag'] for item in pars['hazard']]
    event_occurrence = get_events_occurrence(pars)
    m_min = pars['m_min']
    m_max = pars['m_max']
    lamb = pars['lambda']
    time_periods = pars['time_intervals']
    cpr_min = 1e-6
    d_mag = 0.005
    cpr_v = []
    temp_m_max = round(m_max, 3)
    event_occurrence.m_max = temp_m_max
    m_v = np.arange(m_min, temp_m_max + d_mag / 2, d_mag)
    if what == 'lambda':
        rp_v = []
        for mag in np.nditer(m_v):
            [lambda_mag, rp] = return_period(mag, lamb, event_occurrence.magnitude_distribution)
            rp_v.append(lambda_mag)
        m_color = 'b'
        m_width = 2
        ax.plot(m_v, rp_v, color=m_color, linewidth=m_width)
        ax.set_ylabel('Lambda for magnitude')
    elif what == 'return_period':
        rp_v = []
        for mag in np.nditer(m_v):
            [lambda_mag, rp] = return_period(mag, lamb, event_occurrence.magnitude_distribution)
            rp_v.append(rp)
        m_color = 'b'
        m_width = 2
        ax.plot(m_v, rp_v, color=m_color, linewidth=m_width)
        ax.set_ylabel('Lambda for magnitude')
    elif what == 'probabilities':
        if 'time_intervals' not in pars:
            return
        m_colors = ['k', 'r', 'g', 'b', 'm', 'c', 'y']
        m_width = 2
        for idx, per in enumerate(pars['time_intervals']):
            probabilities = []
            for mag in np.nditer(m_v):
                [lambda_mag, rp] = return_period(mag, lamb, event_occurrence.magnitude_distribution)
                probability = float(event_occurrence.sf(mag, per))
                probabilities.append(probability)
            ax.plot(m_v, probabilities, color=m_colors[idx], linewidth=m_width, label=f"T = {per:.0f} [Y]")
        ax.legend()
        ax.set_ylabel('Probability of not exceedance')
    else:
        raise 'Can not plot {}'.format(what)
    ax.set_xlabel('Magnitude')
    ax.set_title('Area: {}'.format(pars['area_name']))
    ax.grid(which='both', linestyle=':')


def plot_results(configuration):
    fig1, ax = plt.subplots()
    ax.set_yscale("log")
    plot_prob(configuration, 1.0, ax)
    fig1.savefig('annual_probability.png')
    fig1.show()

    fig2, ax = plt.subplots()
    ax.set_yscale("log")
    plot_return_period(configuration, ax)
    fig2.savefig('return_period.png')
    fig2.show()

    # fig4, ax = plt.subplots()
    # ax.set_yscale("log")
    # plot_hazard(configuration, 'lambda', ax)
    # fig4.savefig('lambda.png')
    # fig4.show()

    fig5, ax = plt.subplots()
    # ax.set_yscale("log")
    plot_hazard(configuration, 'probabilities', ax)
    fig5.savefig('probabilities.png')
    fig5.show()
    print("Close all figures to continue")
    plt.show(block=True)


def main():
    params = load_configuration()
    plot_results(params)


if __name__ == "__main__":
    main()
