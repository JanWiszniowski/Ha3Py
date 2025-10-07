r"""
Library for creating the configuration file
-------------------------------------------

The procedure realises Bayesian MAP, Max of Posterior Fiduicial & Gauss (Kijko, 2004)
The algorithm name in the configuration is 'bayesian_fiducial'.

:copyright:
    Jan Wiszniowski <jwisz@igf.edu.pl>,
    Andrzej Kijko <andrzej.kijko@up.ac.za>
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
:version 0.0.1:
    2025-01-01
"""

import sys
import os
import json
import datetime
import numpy as np
from math import floor
from ha3py.GutenbergRichter import GutenbergRichter
# from ha3py.CompoundGutenbergRichter import CompoundGutenbergRichter
from ha3py.utils import HaPyException
from ha3py.constant_values import LN_10_, METHODS
from ha3py.get_events_occurrence import get_events_occurrence


def input_configuration_name():
    if len(sys.argv) >= 2:
        if os.path.isfile(sys.argv[1]):
            return sys.argv[1]
    return None


def load_configuration(required=True):

    if len(sys.argv) >= 2:
        if os.path.isfile(sys.argv[1]):
            try:
                with open(sys.argv[1], "r") as f:
                    return json.load(f)
            except json.decoder.JSONDecodeError as er:
                print(f"The configuration file error: {er}")
        else:
            print(f"The configuration file {sys.argv[1]} does not exist")
    if required:
        print(f"Call: {sys.argv[0]} <params.json>")
        exit(-1)
    return dict()


def save_configuration(configuration, file_name=None):
    if file_name is None:
        file_name = configuration.get('output_configuration', 'Ha3Py.json')
    with open(file_name, 'w') as output_file:
        json.dump(configuration, output_file, indent=4)


def date_years(year, month, day):
    r"""
    Coverts date (year, month, day) into years as float value
    :param year:
    :param month:
    :param day:
    :return: years
    """
    fday = float(day)
    iday = floor(fday)
    fhour = (fday - iday) * 24.0
    ihour = floor(fhour)
    fmin = (fhour - ihour) * 60.0
    imin = floor(fmin)
    fsec = (fmin - imin) * 60.0
    isec = floor(fsec)
    fusec = (fsec - isec) * 1000000.0
    iusec = round(fusec)
    iyear = int(year)
    if iyear >= 1:
        result_date_years = datetime.datetime(int(year), int(month), iday, ihour, imin, isec,
                                              iusec).toordinal() / 365.25
    elif iyear == 0:
        raise HaPyException("Zero year doesn't exist")
    else:
        result_date_years = datetime.datetime(1, int(month), iday, ihour, imin, isec,
                                              iusec).toordinal() / 365.25
        result_date_years += iyear
    return result_date_years
    # int_year = int(year)
    # if int_year % 4 and not int_year % 100:
    #     return float(year) + (YMD_366[int(month) - 1] + float(day) - 1) / 366
    # return float(year) + (YMD_365[int(month) - 1] + float(day) - 1) / 365


def set_if(pars, key, prompt, dtype, ext=''):
    if key in pars:
        return False
    text = input(prompt)
    if dtype == str:
        if ext and text.find(".") == -1:
            text += '.' + ext
        pars[key] = text
    else:
        pars[key] = dtype(text)
    return True


def set_modify(pars, key, prompt, extra, dtype):
    if key in pars:
        wyn = input(prompt.format(extra.format(pars[key])))
        if not wyn:
            return False
        pars[key] = dtype(wyn)
    else:
        pars[key] = dtype(input(prompt.format('')))
    return True


def set_with_warning(pars, key, value):
    if key in pars:
        if pars[key] != value:
            answer = input(f"Parameter '{key}' should be {value} and is {pars[key]} [replace/ignore/cancel] > ")
            if answer[0] == 'i':
                return
            if answer[0] == 'c':
                quit()
    pars[key] = value


def set_if_default(pars, key, prompt, default):
    if key in pars:
        return True
    wyn = input(prompt.format(default))
    if not wyn:
        pars[key] = default
        return True
    wyn = float(wyn)
    if wyn < default:
        pars[key] = default
    else:
        pars[key] = float(wyn)
    return True


def set_modify_default(pars, key, prompt, default):
    return_value = True
    if key in pars:
        default = pars[key]
        return_value = False
    wyn = input(prompt.format(default))
    if  wyn:
        pars[key] = float(wyn)
        return True
    pars[key] = default
    return return_value


def set_uncertainty(pars, key, prompt, default):
    return_value = True
    if key in pars:
        default = 100.0 / np.sqrt(pars[key])
        return_value = False
    wyn = input(f"Define the {prompt} uncertainty in percents (or enter for {default})  [%] > ")
    if  wyn:
        pars[key] = (100 / float(wyn))**2
        return True
    pars[key] = (100 / default)**2
    return return_value


# def set_(pars, key, prompt, default):
#     return_value = True
#     if key in pars:
#         default = pars[key]
#         return_value = False
#     wyn = input(prompt.format(default))
#     if  wyn:
#         pars[key] = float(wyn)
#         return True
#     pars[key] = default
#     return return_value

def read_cat_date(fn):
    dt = fn.readline().split()
    while len(dt) < 3:
        dt = fn.readline().split()
    return date_years(dt[0], dt[1], dt[2])


def read_cat_float(fn):
    while True:
        element = fn.readline()
        try:
            return float(element)
        except ValueError:
            pass


def read_complete_cat(filename):
    catalog = dict()
    with open(filename, "r") as fn:
        catalog['begin'] = read_cat_date(fn)
        catalog['end'] = read_cat_date(fn)
        catalog['time_span'] = catalog['end'] - catalog['begin']
        catalog['m_min'] = read_cat_float(fn)
        magnitude_uncertainty = read_cat_float(fn)
        lines = fn.readlines()
        earthquakes = []
        for line in lines:
            dt = line.split()
            if len(dt) == 5:
                earthquakes.append({'magnitude': float(dt[3]),
                                    'date': date_years(dt[0], dt[1], dt[2]), 'sd': float(dt[4])})
            elif len(dt) == 4:
                earthquakes.append({'magnitude': float(dt[3]), 'date': date_years(dt[0], dt[1], dt[2]),
                                    'sd': magnitude_uncertainty})
            elif len(dt) == 2:
                earthquakes.append({'magnitude': float(dt[0]), 'sd': float(dt[1])})
            elif len(dt) == 1:
                earthquakes.append({'magnitude': float(dt[0]), 'sd': magnitude_uncertainty})
    catalog['earthquakes'] = [m for m in earthquakes if m['magnitude'] >= catalog['m_min']]
    catalog['name'] = filename
    return catalog


def read_paleo_cat(filename):
    catalog = dict()
    m_min = 100.0
    with open(filename, "r") as fn:
        catalog['begin'] = read_cat_date(fn)
        catalog['end'] = read_cat_date(fn)
        catalog['time_span'] = catalog['end'] - catalog['begin']
        lines = fn.readlines()
        earthquakes = []
        prev_date = catalog['begin']
        for line in lines:
            dt = line.split()
            if len(dt) < 3:
                continue
            magnitude = float(dt[2])
            d1 = date_years(dt[0], 1, 1)
            d2 = date_years(dt[1], 1, 1)
            mean_date = (d1 + d2) / 2
            # datediff = d2-d1
            datediff = mean_date - prev_date
            earthquakes.append({'magnitude': magnitude, 'date': mean_date, 'time_span': datediff, 'sd': float(dt[3])})
            prev_date = mean_date
            if m_min > magnitude:
                m_min = magnitude
    catalog['earthquakes'] = earthquakes
    catalog['m_min'] = m_min
    catalog['name'] = filename
    return catalog


def read_hist_cat(filename):
    catalog = dict()
    with open(filename, "r") as fn:
        
        catalog['begin'] = read_cat_date(fn)
        catalog['end'] = read_cat_date(fn)
        catalog['m_min'] = read_cat_float(fn)
        catalog['time_span'] = catalog['end'] - catalog['begin']
        magnitude_uncertainty = read_cat_float(fn)
        lines = fn.readlines()
        earthquakes = []
        prev_date = catalog['begin']
        for line in lines:
            dt = line.split()
            if len(dt)  < 3:
                continue
            event_date = date_years(dt[0], dt[1], dt[2])
            magnitude = float(dt[3])
            time_span = event_date - prev_date
            if len(dt) == 5:
                earthquakes.append({'magnitude': magnitude, 'date': event_date, 'time_span': time_span, 'sd': float(dt[4])})
            elif len(dt) == 4:
                earthquakes.append({'magnitude': magnitude, 'date': event_date, 'time_span': time_span, 'sd': magnitude_uncertainty})
            prev_date = event_date
    catalog['earthquakes'] = [m for m in earthquakes if m['magnitude'] >= catalog['m_min']]
    last_earthquake = catalog['earthquakes'][-1]
    if last_earthquake['date'] < catalog['end']:
        last_earthquake['time_span'] += catalog['end'] - last_earthquake['date']
    catalog['name'] = filename
    return catalog


def init_lambda_beta(configuration, init_beta=True, init_lambda=True, m_max=None):
    r"""
    Initialize :math:`\beta' and :math:`lambda` based on complete catalogs.
    If coefficient 'beta' exists, do not change it

    :param configuration:
    :param init_beta:
    :param init_lambda:
    :param m_max:
    """
    # if 'beta' not configuration:
    #
    if not configuration.get('complete_catalogs',[]):
        return
    if 'beta' in configuration:
        beta = configuration['beta']
    else:
        sum_val = 0.0
        sum_weights = 0.0
        for catalog in configuration['complete_catalogs']:
            m_min = catalog['m_min']
            magnitudes = [earthquake['magnitude'] for earthquake in catalog['earthquakes']]
            no_earthquakes = len(magnitudes)
            beta_for_catalog = 1.0 / (np.mean(magnitudes) - m_min)
            beta_weight = no_earthquakes / beta_for_catalog**2
            sum_val += beta_for_catalog * beta_weight
            sum_weights += beta_weight
        beta = sum_val / sum_weights
        if init_beta:
            configuration['beta'] = sum_val / sum_weights
    if not init_lambda or 'lambda' in configuration:
        return
    sum_val = 0.0
    sum_weights = 0.0
    magnitude_distribution = GutenbergRichter(configuration, beta=beta, m_max=m_max)
    for catalog in configuration['complete_catalogs']:
        time_span = catalog['time_span']
        m_min = catalog['m_min']
        no_earthquakes = len(catalog['earthquakes'])
        lambda_for_catalog = float(no_earthquakes) / time_span / magnitude_distribution.sf(m_min)
        lambda_weight =  1.0 / lambda_for_catalog
        sum_val += lambda_for_catalog * lambda_weight
        sum_weights += lambda_weight
    if init_lambda:
        configuration['lambda'] = sum_val / sum_weights


def define_m_max_assessment(configuration):
    configuration_modified = False
    print('------------------------------------------------------------------')
    print('Assessment of the maximum regional magnitude m_max is based on:')
    print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ')
    for key, method in METHODS.items():
        print(f" {key}: {method}")
    # configuration_modified |= get_if(configuration, 'procedure_id', 'Enter proper number > ', dtype=int)
    configuration_modified |= set_modify(configuration, 'procedure_id', 'Write proper number{} > ',
                                         ' (or enter to confirm {})', dtype=int)
    procedure_id = configuration['procedure_id']
    configuration['m_max_method'] = METHODS[procedure_id]
    if procedure_id == 1:
        set_with_warning(configuration, 'magnitude_distribution', 'Gutenberg-Richter')
        set_with_warning(configuration, 'occurrence_probability', 'Poisson')
        # set_with_warning(configuration, 'delta', 'Gibowicz-Kijko')
        # configuration['m_max_assessment'] = 'solve_delta'
        configuration.pop('delta', None)
        configuration['m_max_assessment'] = 'Gibowicz-Kijko'
        configuration.pop('bayesian_m_max_assessment', None)
        init_lambda_beta(configuration)
    elif procedure_id == 2:
        set_with_warning(configuration, 'magnitude_distribution', 'Compound Gutenberg-Richter')
        set_with_warning(configuration, 'occurrence_probability', 'Poisson-gamma compound')
        # set_with_warning(configuration, 'delta', 'Gibowicz-Kijko')
        # configuration['m_max_assessment'] = 'solve_delta'
        configuration.pop('delta', None)
        configuration['m_max_assessment'] = 'Gibowicz-Kijko'
        configuration.pop('bayesian_m_max_assessment', None)
        init_lambda_beta(configuration)
    elif procedure_id == 3:
        set_with_warning(configuration, 'magnitude_distribution', 'Gutenberg-Richter')
        set_with_warning(configuration, 'occurrence_probability', 'Poisson')
        set_with_warning(configuration, 'delta', 'Kijko-Sellevoll')
        configuration['m_max_assessment'] = 'solve_delta'
        configuration.pop('bayesian_m_max_assessment', None)
        init_lambda_beta(configuration)
    elif procedure_id == 4:
        set_with_warning(configuration, 'magnitude_distribution', 'Compound Gutenberg-Richter')
        set_with_warning(configuration, 'occurrence_probability', 'Poisson-gamma compound')
        set_with_warning(configuration, 'delta', 'Kijko-Sellevoll')
        configuration['m_max_assessment'] = 'solve_delta'
        configuration.pop('bayesian_m_max_assessment', None)
        init_lambda_beta(configuration)
    elif procedure_id == 5:
        set_with_warning(configuration, 'magnitude_distribution', 'Gutenberg-Richter')
        set_with_warning(configuration, 'occurrence_probability', 'Poisson')
        set_with_warning(configuration, 'delta', 'Tate-Pisarenko')
        configuration['m_max_assessment'] = 'solve_delta'
        configuration.pop('bayesian_m_max_assessment', None)
        init_lambda_beta(configuration)
    elif procedure_id == 6:
        set_with_warning(configuration, 'magnitude_distribution', 'Compound Gutenberg-Richter')
        set_with_warning(configuration, 'occurrence_probability', 'Poisson-gamma compound')
        set_with_warning(configuration, 'delta', 'Tate-Pisarenko')
        configuration['m_max_assessment'] = 'solve_delta'
        configuration.pop('bayesian_m_max_assessment', None)
        init_lambda_beta(configuration)
    elif procedure_id == 7:
        set_with_warning(configuration, 'magnitude_distribution', 'Nonparametric gaussian kernel')
        set_with_warning(configuration, 'occurrence_probability', 'Poisson-gamma compound')
        set_with_warning(configuration, 'delta', 'Kijko-Sellevoll')
        configuration['m_max_assessment'] = 'solve_delta'
        configuration.pop('bayesian_m_max_assessment', None)
        no_largest = configuration.get('no_largest_magnitudes')
        if no_largest is None:
            no_largest = input("Number of largest magnitudes (or enter for all) >")
        else:
            no_largest = input(f"Number of largest magnitudes (or enter for {no_largest}) >")
        if no_largest:
            configuration['no_largest_magnitudes'] = int(no_largest)
        init_lambda_beta(configuration, init_beta=False)
    elif procedure_id == 8:
        set_with_warning(configuration, 'magnitude_distribution', 'Compound Gutenberg-Richter')
        set_with_warning(configuration, 'occurrence_probability', 'Poisson-gamma compound')
        set_with_warning(configuration, 'delta', 'Kijko-Sellevoll')
        configuration['m_max_assessment'] = 'solve_delta'
        configuration['bayesian_m_max_assessment'] = 'bayesian_by_shift'
        configuration['bayesian_m_max_estimator'] = 'expected'
        init_lambda_beta(configuration)
    elif procedure_id == 9:
        set_with_warning(configuration, 'magnitude_distribution', 'Compound Gutenberg-Richter')
        set_with_warning(configuration, 'occurrence_probability', 'Poisson-gamma compound')
        set_with_warning(configuration, 'delta', 'Kijko-Sellevoll')
        configuration['m_max_assessment'] = 'solve_delta'
        configuration['bayesian_m_max_assessment'] = 'bayesian_fiducial'
        configuration['bayesian_m_max_estimator'] = 'expected'
        init_lambda_beta(configuration)
    elif procedure_id == 10:
        set_with_warning(configuration, 'magnitude_distribution', 'Compound Gutenberg-Richter')
        set_with_warning(configuration, 'occurrence_probability', 'Poisson-gamma compound')
        set_with_warning(configuration, 'delta', 'Kijko-Sellevoll')
        configuration['m_max_assessment'] = 'solve_delta'
        configuration['bayesian_m_max_assessment'] = 'fixed_value'
        configuration['bayesian_m_max_estimator'] = 'expected'
        init_lambda_beta(configuration)
    else:
        raise HaPyException(f"Wrong id={procedure_id}")
    if configuration.get('magnitude_distribution', '') == 'Compound Gutenberg-Richter':
        set_uncertainty(configuration, 'q_beta', 'Gutenberg-Richter parameter b', 25)
    if configuration.get('occurrence_probability', '') == 'Poisson-gamma compound':
        set_uncertainty(configuration, 'q_lambda', "mean activity rate 'lambda'", 25)

    if configuration.get('bayesian_m_max_assessment', ''):
        configuration_modified |= \
            set_if_default(configuration, 'prior_m_max',
                           'Prior value of maximum possible earthquake magnitude (not less than {}) > ',
                           configuration['m_max_obs'])
        configuration_modified |= set_if(configuration, 'sd_prior_m_max',
                                         'Standard deviation of prior value of m_max > ', dtype=float)
        if configuration['sd_prior_m_max'] > 9.5:
            configuration['sd_prior_m_max'] = 9.5
        if configuration['sd_prior_m_max'] < 0.1:
            configuration['sd_prior_m_max'] = 0.1
    return configuration_modified


def set_global(configuration, catalogue):
    if configuration['begin'] > catalogue['begin']:
        configuration['begin'] = catalogue['begin']
    if configuration['end'] < catalogue['end']:
        configuration['end'] = catalogue['end']
    if configuration.get('m_min', 20) > catalogue['m_min']:
        configuration['m_min'] = catalogue['m_min']
    sd = catalogue.get('sd', 1.0)
    m_max_obs = configuration.get('m_max_obs', -5)
    for mag in catalogue['earthquakes']:
        if m_max_obs < mag['magnitude']:
            m_max_obs = mag['magnitude']
            configuration['m_max_obs'] = mag['magnitude']
            configuration['sd_m_max_obs'] = mag.get('sd', sd)


def define_catalogs(configuration):
    """
    Procedure define catalogues fields if not exists,
    m_min as minimum completeness of catalogues,
    m_max_obs as maximum magnitude in catalogues,
    general date begin as the earliest beginning in catalogues,
    and general date end as the latest end in catalogues

    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict
    :return:
    :rtype:
    """
    configuration_modified = False
    if 'begin' not in configuration:
        configuration['begin'] = np.inf
        configuration['end'] = -np.inf
    # Paleo-catalog
    configuration_modified |= set_if(configuration, 'paleo_catalog',
                                     'Pre-historic (paleo) data file name (or enter if none) > ', dtype=str)
    if configuration['paleo_catalog']:
        if isinstance(configuration['paleo_catalog'], str):
            paleo_cat = read_paleo_cat(configuration['paleo_catalog'])
            configuration['paleo_catalog'] = paleo_cat
        else:
            paleo_cat = configuration['paleo_catalog']
        if not isinstance(paleo_cat, dict):
            raise Exception("Wrong paleo catalog")
        set_global(configuration, paleo_cat)
    else:
        configuration.pop('paleo_catalog', None)
    # Historic catalog
    configuration_modified |= set_if(configuration,
                                     'historic_catalog', 'Historic data file name (or enter if none) > ',
                                     dtype=str)
    if configuration['historic_catalog']:
        if isinstance(configuration['historic_catalog'], str):
            hist_cat = read_hist_cat(configuration['historic_catalog'])
            configuration['historic_catalog'] = hist_cat
        else:
            hist_cat = configuration['historic_catalog']
        if not isinstance(hist_cat, dict):
            raise Exception("Wrong historical catalog")
        set_global(configuration, hist_cat)
    else:
        configuration.pop('historic_catalog', None)
    # Complete catalogs
    if 'complete_catalogs' not in configuration:
        configuration['complete_catalogs'] = []
        configuration_modified = True
        no_complete_catalogues = int(input('Number of complete data catalogue files > '))
        for idx in range(no_complete_catalogues):
            name_comp_file = str(input(f"Name of the #{idx + 1} file with complete data > "))
            comp_cat = read_complete_cat(name_comp_file)
            configuration['complete_catalogs'].append(comp_cat)
            set_global(configuration, comp_cat)
    else:
        if not isinstance(configuration['complete_catalogs'], list):
            raise Exception('HA3: wrong complete data definition')
        for idx, comp_cat in enumerate(configuration['complete_catalogs']):
            if isinstance(comp_cat, str):
                comp_cat = read_complete_cat(comp_cat)
                configuration['complete_catalog'][idx] = comp_cat
        for idx, comp_cat in enumerate(configuration['complete_catalogs']):
            if not isinstance(comp_cat, dict):
                raise Exception(f"Wrong #{idx + 1} complete catalog")
            set_global(configuration, comp_cat)
    configuration_modified |= set_modify(configuration, 'm_max_obs',
                                         'Maximum (EVER!) observed magnitude determination{} > ',
                                         ' not less or equal to {0} (or enter to confirm {0})',
                                         dtype=float)
    configuration_modified |= set_modify(configuration, 'sd_m_max_obs',
                                         'Standard deviation of the maximum observed magnitude{} > ',
                                         ' (or enter to confirm {})', dtype=float)
    configuration['m_max_current'] = configuration['m_max_obs'] + 0.5
    configuration_modified |= set_modify(configuration, 'm_min', 'Minimum value of magnitude{} > ',
                                         ' less or equal to {0} (or enter to confirm {0})', dtype=float)
    if 'time_span' not in configuration:
        time_begin = input(f"Year when time span starts less or equal to {floor(configuration['begin'])}"
                           f" (or enter to confirm it) >")
        if time_begin:
            time_begin = float(time_begin)
            if time_begin < configuration['begin']:
                configuration['begin'] = time_begin
        configuration['time_span'] = configuration['end'] - configuration['begin']
        configuration_modified = True
    return configuration_modified


def define_magnitude_occurrence(configuration):
    configuration_modified = False
    configuration_modified |= set_if(configuration, 'magnitude_distribution',
                                     'Name of the magnitude distribution model > ', dtype=str)
    configuration_modified |= set_if(configuration, 'occurrence_probability',
                                     'Name of the occurrence probability model > ', dtype=str)
    occurrence_probability,_ = get_events_occurrence(configuration)
    for parameter in occurrence_probability.const_coefficients:
        configuration_modified |= set_if(configuration, 'parameter',
                                         f"Value of the {parameter} > ", dtype=float)
    return configuration_modified


def define_configuration(configuration):
    configuration_modified = False
    configuration_modified |= set_if(configuration, 'output_text_file', 'Name of the output text file > ', dtype=str, ext='txt')
    configuration_modified |= set_if(configuration, 'output_configuration',
                                     'Name of the output configuration file > ', dtype=str, ext='json')
    configuration_modified |= set_if(configuration, 'area_name', 'Name of the area > ', dtype=str)
    print('------------------------------------------------------------------')
    print('Define catalogs:')
    print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ')
    define_catalogs(configuration)
    if 'time_intervals' not in configuration:
        configuration_modified = True
        configuration['time_intervals'] = [1.0]
        print('------------------------------------------------------------------')
        print('Define time intervals, for which seismic hazard will be estimated:')
        print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ')
        print('Time interval #1 = 1 year by default. You do not need define thr first time interval')
        configuration['time_intervals'].append(float(input('Time interval #2 > ')))
        configuration['time_intervals'].append(float(input('Time interval #3 > ')))
        configuration['time_intervals'].append(float(input('Time interval #4 > ')))
    # configuration_modified |= get_if_default(params, 'uncertainty_beta',
    #                                          'Uncertainty of Gutenberg-Richter parameter b [%] > ', 1.0)
    # configuration_modified |= get_if_default(params, 'uncertainty_lambda',
    #                                          'Uncertainty of the mean activity rate - lambda [%] > ', 1.0)
    configuration_modified |= define_m_max_assessment(configuration)

    print('------------------------------------------------------------------')
    print('Define remaining parameters:')
    print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ')
    configuration_modified |= set_if(configuration, 'induced_seismicity',
                                     'Is provision for induced seismicity required (yes/no)? > ', dtype=str)
    configuration['induced_seismicity'] = configuration['induced_seismicity'].lower()
    if configuration['induced_seismicity'] != 'yes':
        configuration['induced_seismicity'] = 'no'
    if configuration['induced_seismicity'].lower() == 'yes':
        configuration_modified |= set_if(configuration, 'induced_seismicity_coefficient',
                                         'Multiplicative factor of activity rate (lambda) > ', dtype=float)
    else:
        configuration['induced_seismicity_coefficient'] = 1.0
    if 'prior_beta' not in configuration:
        configuration_modified |= set_if(configuration, 'prior_b',
                                         'Prior value of the Gutenberg-Richter parameter b (or enter if not defined) > ',
                                         dtype=str)
        if not configuration.get('prior_b', ''):
            configuration.pop('prior_beta', None)
            configuration.pop('sd_prior_beta', None)
            configuration.pop('sd_prior_b', None)
        else:
            configuration_modified |= set_if(configuration, 'sd_prior_b',
                                             'Standard deviation of the prior value of b > ', dtype=float)
            configuration['prior_beta'] = float(configuration['prior_b']) * LN_10_
            configuration['sd_prior_beta'] = configuration['sd_prior_b'] * LN_10_
    configuration_modified |= set_if(configuration, 'likelihood_optimization_method',
                                     'Choose optimization method [see scipy.optimize] (or enter for default)> ',
                                     dtype=str)
    if not configuration['likelihood_optimization_method']:
        configuration['likelihood_optimization_method'] = None
    if configuration_modified:
        configuration['created_on'] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return configuration_modified

def main():
    print("==============================================================")
    configuration = load_configuration(required=False)
    if configuration:
        if 'area_name' in configuration:
            print(f"Loading configuration file for '{configuration['area_name']}'")
        else:
            print(f"Loading configuration file for undefined area")
        print("==============================================================")
    else:
        print("The configuration will be created from the beginning")
        print("==============================================================")
    if define_configuration(configuration):
        if len(sys.argv) >= 3:
            save_configuration(configuration, file_name=sys.argv[2])
        else:
            save_configuration(configuration)


if __name__ == "__main__":
    main()
