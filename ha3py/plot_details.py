import sys
import numpy as np
# import matplotlib
import matplotlib.pyplot as plt
from ha3py.configuration import load_configuration, set_if_default
from ha3py.Delta import get_delta
from ha3py.get_events_occurrence import get_events_occurrence
from ha3py.get_magnitude_distribution import get_magnitude_distribution


# from ha3py.print_results import print_results
# from ha3py.ln_likelihood import ln_likelihood


def plot_catalog(catalog, bin_width, title, axis):
    magnitudes = np.array([event['magnitude'] for event in catalog['earthquakes']])
    m_min = catalog['m_min']
    bins = np.arange(m_min - bin_width / 2.0, np.max(magnitudes) + bin_width, bin_width)
    axis.hist(magnitudes, bins=bins, log=True)
    axis.set_title(f"{title}, {catalog['begin']:6.2f}-{catalog['end']:6.2f}")
    # print(f"{title}, {catalog['begin']:6.2f}-{catalog['end']:6.2f}")
    y_lim = axis.get_ylim()
    # axis.semilogy((m_min, m_min), y_lim, 'r')
    if y_lim[1] > 10:
        return
    yticklabels = axis.get_ymajorticklabels()
    for label in yticklabels:
        y = label.get_position()[1]
        if y == 1.0:
            label.set_text(f"{1}")
    axis.set_yticklabels(yticklabels)
    yticklabels = axis.get_yminorticklabels()
    for label in yticklabels:
        y = label.get_position()[1]
        if 1 < y < 11 and label.get_text():
            # print(label)
            label.set_text(f"{y:.0f}")
    axis.set_yticklabels(yticklabels, minor=True)


def catalogs_histogram(configuration):
    set_if_default(configuration, 'plot_magnitudes_bin_width', "Histogram's bins width >", 0.2)
    bin_width = configuration['plot_magnitudes_bin_width']
    no_catalogs = 0
    p_phs = configuration.get('paleo_catalog')
    p_his = configuration.get('historic_catalog')
    p_comp = configuration.get('complete_catalogs')
    if p_phs:
        no_catalogs += 1
    if p_his:
        no_catalogs += 1
    if p_comp:
        no_catalogs += len(p_comp)
    figure, axes = plt.subplots(no_catalogs, 1, sharex=True, squeeze=False)
    axis_index = 0
    if p_phs:
        plot_catalog(p_phs, bin_width, 'Paleo-catalog', axes[axis_index][0])
        axis_index += 1
    if p_his:
        plot_catalog(p_his, bin_width, 'Historical catalog', axes[axis_index][0])
        axis_index += 1
    if p_comp:
        for idx, comp in enumerate(p_comp):
            plot_catalog(comp, bin_width, f'#{idx + 1} complete catalog', axes[axis_index][0])
            axis_index += 1
    axes[-1][0].set_xlabel('Magnitude')


def plot_delta(configuration):
    r"""
    It plots :math:`\Delta + m_{max}^{obs}` as a function of :math:`m_{max}`.
    For notice the solution :math:`m_{max}` is plotted as dotted line

    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict

    """
    number_of_points = 100
    event_occurrence = get_events_occurrence(configuration)
    if 'm_min_ref' in configuration:
        event_occurrence.m_min = configuration['m_min_ref']
    lambda_ref = event_occurrence.d_mean()
    delta = get_delta(configuration)
    m_max_obs = configuration['m_max_obs']
    print(f"Maximum observed magnitude {m_max_obs}")
    m_max_abs = m_max_obs + 1
    deltas = np.zeros(number_of_points, dtype=float)
    m_maxes = np.linspace(m_max_obs + 0.01, m_max_obs + 1, number_of_points)
    for idx, m_max in enumerate(m_maxes):
        delta.magnitude_distribution.m_max = m_max
        deltas[idx] = delta(time=configuration['time_span'], annual_lambda=lambda_ref) + m_max_obs
    figure, axis = plt.subplots()
    axis.plot(m_maxes, deltas, color='k', linewidth=2)
    axis.plot([m_max_obs, m_max_abs], [m_max_obs, m_max_abs], 'k--')
    axis.set_ylabel(r"$m_{max}^{obs} + \Delta$")
    axis.set_xlabel(r"$m_{max}$")
    axis.set_title(configuration['delta'])
    # axis.set_ylim(np.min(deltas), np.max(deltas))
    # print_results(configuration)


def plot_non_occurrence(configuration):
    number_of_points = 100
    event_occurrence = get_events_occurrence(configuration)
    magnitudes = np.linspace(event_occurrence.m_min - 0.1, event_occurrence.m_max + 0.1, number_of_points)
    figure, axis = plt.subplots(2, len(configuration['time_intervals']))
    for idx, time_interval in enumerate(configuration['time_intervals']):
        nom_pdf = [event_occurrence.pdf(m, time_interval) for m in magnitudes]
        nom_cdf = [event_occurrence.cdf(m, time_interval) for m in magnitudes]
        axis[0][idx].plot(magnitudes, nom_pdf, color='k', linewidth=2)
        axis[1][idx].plot(magnitudes, nom_cdf, color='k', linewidth=2)
        axis[1][idx].set_xlabel("magnitide")
        axis[0][idx].set_title(f"{configuration['occurrence_probability']},\nnon-occurrence,\nT = {time_interval}",
                               loc='center')
    axis[0][0].set_ylabel("PDF")
    axis[1][0].set_ylabel("CDF")


def plot_magnitude_distribution(configuration):
    r"""
    It plots magnitude distributions CDF and PDF according the configuration.

    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict
    """
    magnitude_distribution = get_magnitude_distribution(configuration)
    magnitudes = np.linspace(magnitude_distribution.m_min - 0.1, magnitude_distribution.m_max + 0.1, 250)
    m_pdf = np.zeros(250, dtype=float)
    m_cdf = np.zeros(250, dtype=float)
    for idx, m in enumerate(magnitudes):
        m_pdf[idx] = magnitude_distribution.pdf(m)
        m_cdf[idx] = magnitude_distribution.cdf(m)
    figure, axis = plt.subplots(2, 1)
    axis[0].plot(magnitudes, m_pdf, color='k', linewidth=2)
    axis[1].plot(magnitudes, m_cdf, color='k', linewidth=2)
    axis[0].set_ylabel("PDF")
    axis[1].set_ylabel("CDF")
    axis[1].set_xlabel("magnitide")
    axis[0].set_title(configuration.get('magnitude_distribution', 'Nonparametric Gaussian kernel'))


def log_magnitude_distribution(configuration):
    r"""
    It plots logarithm of magnitude distributions CDF and PDF according the configuration.

    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict
    """
    magnitude_distribution = get_magnitude_distribution(configuration)
    magnitudes = np.linspace(magnitude_distribution.m_min - 0.1, magnitude_distribution.m_max + 0.1, 250)
    m_pdf = np.zeros(250, dtype=float)
    m_cdf = np.zeros(250, dtype=float)
    for idx, m in enumerate(magnitudes):
        m_pdf[idx] = magnitude_distribution.pdf(m)
        m_cdf[idx] = magnitude_distribution.cdf(m)
    figure, axis = plt.subplots(2, 1)
    axis[0].semilogy(magnitudes, m_pdf, color='k', linewidth=2)
    axis[1].semilogy(magnitudes, m_cdf, color='k', linewidth=2)
    axis[0].set_ylabel("pdf")
    axis[1].set_ylabel("pdf")
    axis[1].set_xlabel("magnitide")
    axis[0].set_title(configuration.get('magnitude_distribution', 'Nonparametric Gaussian kernel'))


def plot_likelihood(configuration):
    pass
    # points = 20
    # val_range = 4
    # print(f"Beta: {configuration['beta']}, lambda {configuration['lambda']}")
    # event_occurrence = get_events_occurrence(configuration)
    # x_v = event_occurrence.coefficients[:-1]
    # print(f"Beta: {x_v[1]}, lambda {x_v[0]}")
    # lambda_values = np.linspace(x_v[0] / val_range, x_v[0] * val_range, points)
    # beta_values = np.linspace(x_v[1] / val_range, x_v[1] * val_range, points)
    # ll = np.zeros((points, points), dtype=float)
    # min_ll = 1e20
    # beta_min = 0
    # lambda_min = 0
    # for beta_index, beta_value in enumerate(beta_values):
    #     x_v[1] = beta_value
    #     print(f"{beta_index}: {beta_value}")
    #     for lambda_index, lambda_value in enumerate(lambda_values):
    #         x_v[0] = lambda_value
    #         val = ln_likelihood(x_v, configuration)
    #         if val < min_ll:
    #             min_ll = val
    #             beta_min = beta_value
    #             lambda_min = lambda_value
    #         ll[beta_index][lambda_index] = val
    # figure, axis = plt.subplots()
    # pcm = axis.pcolor(lambda_values, beta_values, ll, cmap='hsv')
    # figure.colorbar(pcm, ax=axis)
    # print(f"Min. beta: {beta_min}, min. lambda {lambda_min}")


WhatToPlot = {
    'delta': plot_delta,
    'magnitude_distribution': plot_magnitude_distribution,
    'log_magnitude_distribution': log_magnitude_distribution,
    'catalogs_histogram': catalogs_histogram,
    'likelihood': plot_likelihood,
    'non_occurrence': plot_non_occurrence
}


def plot_figure(figure, configuration):
    plotter = WhatToPlot.get(figure)
    if plotter is not None:
        plotter(configuration)
    print("Close the figure to continue")
    plt.show(block=True)


def main():
    if len(sys.argv) < 3:
        print(f"Call: {sys.argv[0]} <configuration.json> <what_to_plot>")
        print(f"  where <what_to_plot> is a list (separated by commas, without spaces) of:")
        for plot_type in WhatToPlot.keys():
            print(f"  - {plot_type}")
        exit(-1)
    configuration = load_configuration()
    figures = sys.argv[2].split(',')
    for figure in figures:
        plotter = WhatToPlot.get(figure)
        if plotter is not None:
            plotter(configuration)
    plt.show(block=True)


if __name__ == "__main__":
    main()
