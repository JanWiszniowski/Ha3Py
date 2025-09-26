"""
Ha3Py
(c) Andrzej Kijko, Jan Wiszniowski
ver. 2024-01
"""

import numpy as np
from numpy.linalg import inv
from ha3py.return_period import return_period
from ha3py.get_events_occurrence import get_events_occurrence
from ha3py.uncertainty import get_covariance
from ha3py.configuration import load_configuration, save_configuration
from ha3py.utils import print_separation_double_line, print_separation_single_line
from ha3py.constant_values import LN_10_


def print_share(prompt, share, parameter_names):
    prompt += ':'
    for idx, name in enumerate(parameter_names):
        print(f'{prompt:<23} {name:<7}= {share[idx][idx]:4.1f}%')
        prompt = ' '


def compute_hazard(configuration):
    event_occurrence = get_events_occurrence(configuration)
    lamb = configuration['lambda']
    m_min = configuration['m_min']
    m_max = configuration['m_max']
    time_periods = configuration['time_intervals']
    d_mag = 0.1
    # =========================================================================
    mag_start = round(m_min, 1)
    # mag_end = round(m_max, 1)
    mag_end = m_max
    haz_a = []
    if mag_end > mag_start:
        for mag in np.nditer(np.arange(mag_start, mag_end, d_mag)):  # MAGNITUDE LOOP STARTS HERE >>>>>>>>
            [lambda_mag, rp] = return_period(mag, lamb, event_occurrence.magnitude_distribution)
            probabilities = [float(event_occurrence.sf(mag, tp)) for tp in time_periods]
            haz_a.append({'mag': float(mag), 'lambda_mag': float(lambda_mag), 'return_period': float(rp), 'probabilities': probabilities})
        # END OF MAGNITUDE LOOP <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    return haz_a
    # END OF FUNCTION compute_hazard ==========================================
    # =========================================================================


def print_results(configuration):
    print_separation_double_line()
    print('Final results')
    print_separation_single_line()
    print('Area: {}'.format(configuration['area_name']))
    print('Created on {}'.format(configuration['created_on']))
    print('Computed on {}'.format(configuration['computation_time']))
    print_separation_single_line()
    print('Occurrence probability: '.format(configuration['occurrence_probability']))
    print('Magnitude distribution: '.format(configuration['magnitude_distribution']))
    print('M_max assessment method: '.format(configuration['m_max_method']))
    print_separation_single_line()
    cov = configuration['cov_beta_lambda']
    events_distribution_names = configuration['coefficient_names']

    for name in events_distribution_names:
        if name == 'beta':
            print(f"{name:<10s} = {configuration.get(name, -999):7.3f} +/- {configuration.get('sd_' + name):.3f},"
                  f" (b = {configuration.get(name, -999)/LN_10_:.2f} "
                  f"+/- {configuration.get('sd_' + name)/LN_10_:.2f})")
        elif name == 'lambda':
            print(f"{name:<10s} = {configuration.get(name, -999):7.3f} +/- {configuration.get('sd_' + name):.3f},"
                  f" (for m_min = {configuration['m_min']:.2f})")
            if 'induced_seismicity' in configuration and configuration['induced_seismicity'] == 'yes':
                print(f"    LAMBDA(IS) = {configuration['lambda_is']:.3f} +/-{configuration['sd_lambda_is']:.3f}")
        else:
            print(f"{name:<10s} = {configuration.get(name, -999):7.3f} +/- {configuration.get('sd_' + name):.3f}")
    # "m_max_method": "Kijko-Sellevoll-Bayes (Kijko, 2004)",
    # "magnitude_distribution": "Compound Gutenberg-Richter",
    # "q_beta": 16,
    # "occurrence_probability": "Poisson-gamma compound"
    #  print(' It was calculated according to {} method'.format(configuration.get('m_max_method',
    #                                                                            'manually defined')))
    print_separation_single_line()
    first_line = True
    for cov_row in cov:
        if first_line:
            cov_str ='COV = [ '
            first_line = False
        else:
            cov_str = '      [ '
        for value in cov_row:
            cov_str += f'{value:8.3f} '
        cov_str += ']'
        print(cov_str)
    print('')
    for idx1 in range(len(events_distribution_names)-1):
        for idx2 in range(idx1+1, len(events_distribution_names)-1):
            print("Corr({},{}) = {:.3f}".format(events_distribution_names[idx1], events_distribution_names[idx2],
                                                cov[idx1][idx2] / np.sqrt(cov[idx1][idx1] * cov[idx2][idx2])))
    print('')


def print_hazard(configuration):
    s = '| Mag  |Lambda(sf)|   RP    |'
    for rp in configuration['time_intervals']:
        s += ' pr. T={:5.0f}|'.format(rp)
    line_eq = '=' * len(s)
    print(line_eq)
    spaces_before = (len(s)-len('SEISMIC HAZARD')-2)//2
    title = '|' + ' ' * spaces_before + 'SEISMIC HAZARD'
    title = title + ' ' * (len(s) - len(title) - 1) + '|'
    print(title)
    print(line_eq)
    print(s)
    s = '|------|----------|---------|'
    for rp in configuration['time_intervals']:
        s += '------------|'.format(rp)
    print(s)
    for hz in configuration['hazard']:
        s = '| {:4.2f} | {:6.2e} | {:7.1f} |'.format(hz['mag'], hz['lambda_mag'], hz['return_period'])
        for pr in hz['probabilities']:
            s += '  {:8.6f}  |'.format(pr)
        print(s)
    print(line_eq)


def print_percentage_share(configuration):
    p_phs = configuration.get('paleo_catalog')
    p_his = configuration.get('historic_catalog')
    p_comp = configuration.get('complete_catalogs')
    events_distribution = get_events_occurrence(configuration)
    names = events_distribution.coefficient_names[:-1]
    # x_v = events_distribution._coefficient_values()
    temp_pars = configuration.copy()
    temp_pars.pop('paleo_catalog', None)
    temp_pars.pop('historic_catalog', None)
    temp_pars.pop('complete_catalogs', None)
    # =========================================================================
    # Calculations of info provided by PALEO part of catalogue ================
    if p_phs:
        temp_pars['paleo_catalog'] = p_phs
        var_cov_a = get_covariance(temp_pars)
        temp_pars.pop('paleo_catalog', None)
        info_phs_a = abs(inv(var_cov_a))
    else:
        info_phs_a = np.zeros((2, 2))
    # End of info assessment provided by PALEO part of catalogue ==============
    # =========================================================================
    # Calculations of info provided by HISTORIC part of catalogue =============
    if p_his:
        temp_pars['historic_catalog'] = p_his
        var_cov_a = get_covariance(temp_pars)
        temp_pars.pop('historic_catalog', None)
        info_his_a = abs(inv(var_cov_a))
    else:
        info_his_a = np.zeros((2, 2))
    # End of info assessment provided by HISTORIC catalogue ===================
    # =========================================================================
    # Calculations of info provided by COMPLETE part of catalogue =============
    if p_comp:
        temp_pars['complete_catalogs'] = p_comp
        var_cov_a = get_covariance(temp_pars)
        temp_pars.pop('complete_catalogs', None)
        info_comp_a = abs(inv(var_cov_a))
    else:
        info_comp_a = np.zeros((2, 2))
    # End of info assessment provided by COMPLETE catalogue ===================
    # =========================================================================
    # Info by WHOLE & percent of info provided by EACH catalogue ==============
    info_whole_a = info_phs_a + info_his_a + info_comp_a  # WHOLE CATALOGUE
    info_phs_v = 100.0 * info_phs_a / info_whole_a
    info_his_v = 100.0 * info_his_a / info_whole_a
    info_comp_v = 100.0 * info_comp_a / info_whole_a
    # =========================================================================
    # DISPLAY INFO PROVIDED BY EACH PART OF CATALOGUE =========================
    print_separation_double_line()
    print('Information provided by each part of catalogue (in per-cent)')
    print_separation_single_line()
    if p_phs:
        print_share('Prehistoric catalogue', info_phs_v, names)
    if p_his:
        print_share('Historic catalogue', info_his_v, names)
    if p_comp:
        print_share('Complete catalogue(s)', info_comp_v, names)

    # END OF FUNCTION : print_percentage_share ================================
    # =========================================================================


def main():
    configuration = load_configuration()
    print_results(configuration)
    configuration['hazard'] = compute_hazard(configuration)
    print_hazard(configuration)
    save_configuration(configuration)


if __name__ == "__main__":
    main()

