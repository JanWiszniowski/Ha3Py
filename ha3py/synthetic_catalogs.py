"""
The synthetic catalogs simulates the paleo-, historical, and complete catalogues.
The package is the Python conversion
of the MATLAB MC_data_generation package by Andrzej Kijko.

..

    :copyright:
        Jan Wiszniowski <jwisz@igf.edu.pl>,
        Andrzej Kijko <andrzej.kijko@up.ac.za>
    :license:
        GNU Lesser General Public License, Version 3
        (https://www.gnu.org/copyleft/lesser.html)
    :version 0.0.3:
        2026-03-03


"""

from ha3py.utils import HaPyException
from ha3py.configuration import load_configuration, save_configuration
from ha3py.constant_values import EPS
from ha3py.get_magnitude_distribution import get_magnitude_distribution
import numpy as np


def test_lambda_beta(catalog_configuration, printing=True):
    r"""
    The function estimated and prints :math:`\beta` and :math:`\lambda` parameters
    of the catalogue assessed by the simplest method:

    .. math::
        \beta =\frac{1}{n}\sum_{i=1}^{n}m_i -m_{min},

    and

    .. math::
        \lambda = \frac{n}{T},

    where :math:`T` is the catalogue time span, :math:`n` is the number of events,
    and :math:`m_{min}` the minimum magnitude of completeness.

    :param catalog_configuration: The catalog definition, part of the dictionary of all Ha3Py parameters
    :type catalog_configuration: dict
    :return: None
    """
    m_min = catalog_configuration['m_min']
    magnitudes = [earthquake['magnitude'] for earthquake in catalog_configuration['earthquakes']]
    time_span = catalog_configuration['time_span']
    no_earthquakes = len(magnitudes)
    beta_catalog = 1.0 / (np.mean(magnitudes) - m_min)
    # magnitude_distribution = GutenbergRichter(configuration, beta=beta_for_catalog)
    lambda_catalog = float(no_earthquakes) / time_span
    if printing:
        print(f"Catalog {catalog_configuration['name']}: beta={beta_catalog}, lambda={lambda_catalog} for m_min={m_min}")
    return beta_catalog, lambda_catalog


def f_rnd_app_mag(mag, sd_mag):
    # FUNCTION PRODUCES APPARENT VALUE OF EARTHQUAKE MAGNITUDE BY ADDING
    # GAUSSIAN (NORMAL) ERROR WITH ZERO MEAN AND STANDARD DEVIATION SD_MAG
    # TO THE GIVEN EARTHQUAKE MAGNITUDE MAG
    #
    #
    # WRITTEN           : 23 FEB 1999
    # LAST TIME REVISED : 10 DEC 2000
    #
    #
    # INPUT PARAMETERS:
    #
    #   mag    = EARTHQUAKE MAGNITUDE
    #   sd_mag = STANDARD DEVIATION OF GAUSSIAN ERROR OF EARTHQUAKE
    #            MAGNITUDE DETERMINATION
    #
    # ======================================================================
    sd_mag3 = 3.0 * sd_mag
    d_mag = np.random.normal() * sd_mag  # NORMAL (GAUSSIAN) ERROR OF MAGNITUDE
    if abs(d_mag) > sd_mag3:
        if d_mag > 0:
            d_mag = sd_mag3
        if d_mag < 0:
            d_mag = -sd_mag3
    app_mag = mag + d_mag  # APPARENT MAGNITUDE
    return 0.01 * round(100.0 * app_mag)
    # END OF rnd_app_mag_f =================================================


def f_rnd_mag_grb(configuration):
    # FUNCTION GENERATES RANDOM VALUE OF EARTHQUAKE MAGNITUDE
    # FOLLOWING THE GUTENBERG-RICHTER-BAYESIAN DISTRIBUTION
    #
    # WRITTEN           : 09 DEC 2000
    # LAST TIME REVISED : 10 DEC 2000
    #
    #
    # INPUT PARAMETERS: p_grb_v = vector of parameters of Gutenberg-Richter-
    #                             BAYESIAN magnitude distribution, where
    #
    #   p_grb_v(1,1) = beta parameter of Gutenberg-Richter
    #   p_grb_v(2,1) = (beta / sd_beta)^2
    #   p_grb_v(3,1) = m_min = level of completeness
    #   p_grb_v(4,1) = m_max
    #
    # ======================================================================

    beta = configuration['beta']  # beta parameter of Gutenberg-Richter
    q_beta = configuration.get('q_beta', 10000)  # (beta / sd_beta)^2
    m_min = configuration['m_min']  # level of completeness
    m_max = configuration['m_max']  # m_max
    if q_beta < 1.0:
        q_beta = 1.0
    if q_beta > 10000:
        q_beta = 10000
    tmp = q_beta / (q_beta + beta * (m_max - m_min))
    tmp = 1 - tmp ** q_beta
    c = 1 / tmp  # NORMALIZING COEFFICIENT OF G-R-B CDF
    r_n = np.random.rand()  # UNIFORMLY DISTRIBUTED RANDOM NUMBER
    if r_n < EPS:
        r_n = EPS
    if (1.0 - r_n) < EPS:
        r_n = 1.0 - EPS
    tmp = q_beta / ((1.0 - r_n / c) ** (1 / q_beta)) - q_beta
    mag = m_min + tmp / beta
    if mag < m_min:
        mag = m_min + EPS
    if mag > m_max:
        mag = m_max - EPS
    # return mag
    return 0.01 * round(100 * mag)
    # END OF rd_mag_grb_f =================================================


def f_rnd_time_int(lamb):
    # FUNCTION GENERATES RANDOM VALUE OF TIME INTERVAL BETWEEN SUCCESSIVE EQ-s
    # DISTRIBUTED ACCORDING TO EXPONENTIAL DISTRIBUTION WITH PARAMETER LAMBDA
    #
    # WRITTEN           : 23 FEB 1999
    # LAST TIME REVISED : 07 OCT 2000
    #
    #
    # INPUT PARAMETER: LAMBDA = parameter of the EXPONENTIAL distribution
    #
    # ======================================================================
    if lamb <= 0:
        raise HaPyException(' INPUT parameter LAMBDA must be positive')
    eps3 = 3.0 * EPS
    r_n = np.random.rand()
    if r_n < eps3:
        r_n = eps3
    if (1.0 - r_n) < eps3:
        r_n = 1.0 - eps3
    y = -np.log(1.0 - r_n) / lamb  # RANDOM TIME INTERVAL
    return y
    # END OF f_rnd_time_int ================================================


def f_rnd_poisson(lamb):
    r = 0
    p = 0
    while True:
        p = p - np.log(np.random.rand())
        if p < lamb:
            r += 1
        else:
            break
    return r


def get_times(catalog_configuration):
    """
    The support function returns begin, end time of the catalogue,
    and the time span of the catalogue in years.
    At least two of these values must be defined in the configuration.

    :param catalog_configuration: The catalog definition, part of the dictionary of all Ha3Py parameters
    :type catalog_configuration: dict
    :return: The beginning, end, and time span of teh catalog definition as float years
    :rtype: tuple
    """
    begin_time = catalog_configuration.get('begin')
    end_time = catalog_configuration.get('end')
    time_span = catalog_configuration.get('time_span')
    if begin_time is None:
        begin_time = end_time - time_span
    if end_time is None:
        end_time = begin_time + time_span
    if time_span is None:
        time_span = end_time - begin_time
    return begin_time, end_time, time_span


def prehistorical_simulation(catalog_configuration, configuration):
    # PRE-HISTORIC DATA COMPUTATIONS --------------------------------------
    name = catalog_configuration.get('name', 'paleo_synthetic')
    nr_phs_eq = 0
    m_max_obs = 0.0

    begin_time, end_time, time_span = get_times(catalog_configuration)
    lamb = configuration['lambda']
    m_max = configuration['m_max']
    m_min_current = catalog_configuration.get('m_min', m_max - 1.5)
    time_interval = catalog_configuration.get('time_interval', 100.0)
    time_uncertainty = catalog_configuration.get('time_uncertainty', 50.0)
    sd_mag = catalog_configuration.get('sd', 0.5)
    time_current = begin_time
    last_time_m_max_obs_current = time_current
    earthquakes = []
    while time_current <= end_time:
        lt = lamb * time_interval
        nr_eq_current = f_rnd_poisson(lt)  # CURRENT, RANDOM NR OF EQ-s WITHIN CURRENT time_interval
        if nr_eq_current > 0:
            m_max_obs_current = 0.0
            for i_eq in range(nr_eq_current):
                mag = f_rnd_mag_grb(
                    configuration)  # EARTHQUAKE MAGNITUDE, FROM G-R OR G-R-B DISTRIBUTIONS, STARTING FROM Mmin TOTAL
                if mag < m_min_current:
                    continue
                if mag >= m_max_obs_current:
                    m_max_obs_current = mag  # MAX OBS MAG WITHIN CURRENT TIME INTERVAL time_interval
            m_max_obs_current = f_rnd_app_mag(m_max_obs_current, sd_mag)
            if m_max_obs_current:
                time_m_max_obs_current = time_current + np.random.rand() * time_interval
                # year1 = time_m_max_obs_current - time_uncertainty
                # year2 = time_m_max_obs_current + time_uncertainty
                # year1 = round(year1)  # LOWER BOUND OF EQ-e OCCURRENCE
                # year2 = round(year2)  # UPPER BOUND OF EQ-e OCCURRENCE
                earthquakes.append({"magnitude": m_max_obs_current, "date": time_m_max_obs_current,
                                    "time_span": time_m_max_obs_current - last_time_m_max_obs_current, "sd": sd_mag})
                if m_max_obs_current > m_max_obs:
                    m_max_obs = m_max_obs_current
                nr_phs_eq += 1
            time_current = time_current + time_interval
    print(f"NUMBER of PRE-HISTORIC EQ-s (Mag >= {m_min_current}) = {nr_phs_eq}")
    print(f"PRE-HISTORIC EQ-s are selected from time interval = {time_interval} [Y]")
    print(f"Maximum OBSERVED magnitude (APPARENT) = {m_max_obs}")
    return {'begin': begin_time, 'end': end_time, 'time_span': time_span, 'm_min': m_min_current,
            'sd': time_uncertainty, 'name': name, 'earthquakes': earthquakes}, m_max_obs, m_min_current


def historical_simulation(catalog_configuration, configuration):
    # HISTORIC DATA COMPUTATIONS --------------------------------------
    name = catalog_configuration.get('name', 'historic_synthetic')
    nr_his_eq = 0
    m_max_obs = 0.0
    begin_time, end_time, time_span = get_times(catalog_configuration)
    lamb = configuration['lambda']
    m_max = configuration['m_max']
    m_min_current = catalog_configuration.get('m_min', m_max - 2.0)
    time_interval = catalog_configuration.get('time_interval', 10.0)
    sd_mag = catalog_configuration.get('sd', 0.5)
    time_current = begin_time
    last_time_m_max_obs_current = time_current
    earthquakes = []
    while time_current <= end_time:
        lt = lamb * time_interval
        nr_eq_current = f_rnd_poisson(lt)  # CURRENT, RANDOM NR OF EQ-s WITHIN CURRENT time_interval
        if nr_eq_current > 0:
            m_max_obs_current = 0.0
            for i_eq in range(nr_eq_current):
                mag = f_rnd_mag_grb(
                    configuration)  # EARTHQUAKE MAGNITUDE, FROM G-R OR G-R-B DISTRIBUTIONS, STARTING FROM Mmin TOTAL
                if mag < m_min_current:
                    continue
                if mag >= m_max_obs_current:
                    m_max_obs_current = mag  # MAX OBS MAG WITHIN CURRENT TIME INTERVAL time_interval
            m_max_obs_current = f_rnd_app_mag(m_max_obs_current, sd_mag)
            if m_max_obs_current:
                time_m_max_obs_current = time_current + np.random.rand() * time_interval
                earthquakes.append({"magnitude": m_max_obs_current, "date": time_m_max_obs_current,
                                    "time_span": time_m_max_obs_current - last_time_m_max_obs_current, "sd": sd_mag})
                if m_max_obs_current > m_max_obs:
                    m_max_obs = m_max_obs_current
                nr_his_eq += 1
            time_current = time_current + time_interval
    print(f"NUMBER of HISTORIC EQ-s (Mag >= {m_min_current}) = {nr_his_eq}")
    print(f"HISTORIC EQ-s are selected from time interval = {time_interval} [Y]")
    print(f"Maximum OBSERVED magnitude (APPARENT) = {m_max_obs}")
    return {'begin': begin_time, 'end': end_time, 'time_span': time_span, 'm_min': m_min_current,
            'sd': sd_mag, 'name': name, 'earthquakes': earthquakes}, m_max_obs, m_min_current


def complete_simulation(catalog_configuration, configuration, idx):
    # COMPLETE DATA ========================================================
    # WRITE "HEAD" OF COMPLETE DATA OUTPUT FILE -------------------------
    #  fprintf(id_data_file,'%5d %3d %3d\n',year_begin,month_begin,day_begin);
    #  fprintf(id_data_file,'%5d %3d %3d\n',year_end,month_end,day_end);
    #  fprintf(id_data_file,'%5.2f\n',m_min_current);
    #  fprintf(id_data_file,'%5.2f\n',sd_mag);
    # END OF WRITING OF "HEAD" OF COMPLETE OUTPUT DATA FILE -------------

    # COMPLETE DATA COMPUTATIONS ----------------------------------------
    name = catalog_configuration.get('name', f'complete#{idx}_synthetic')
    time = 0.0
    nr_eq_total = 0
    nr_eq_current = 0
    m_max_obs = 0.0
    mean_mag = 0.0
    m_min_current = catalog_configuration['m_min']
    m_min = catalog_configuration['m_min']
    lamb = configuration['lambda']
    # config_m_min = configuration['m_min']
    sd_mag = catalog_configuration.get('sd', 0.5)
    begin_time, end_time, time_span = get_times(catalog_configuration)
    earthquakes = []
    if sd_mag > (m_min_current - m_min) / 3.0:
        sd_mag = (m_min_current - m_min) / 3.0
    # if m_min_current < config_m_min + sd_mag:
    #     configuration['m_min'] = m_min_current - sd_mag
    #     magnitude_distribution = get_magnitude_distribution(configuration)
    #     catalog_lambda = configuration['lambda'] / magnitude_distribution.sf(config_m_min)
    # else:
    #     catalog_lambda = configuration['lambda']
    while time < time_span:
        dt = f_rnd_time_int(lamb)  # RANDOM TIME INTERVAL [Y] CORRESPONDING TO LAMBDA_TOTAL
        time = time + dt  # CURRENT TIME
        nr_eq_total = nr_eq_total + 1  # NR OF EQ-s STARTING FROM m_min_total
        mag = f_rnd_mag_grb(
            configuration)  # EARTHQUAKE MAGNITUDE, FROM G - R, OR G-R-B DISTRIBUTIONS, STARTING FROM Mmin TOTAL
        mag = f_rnd_app_mag(mag, sd_mag)  # APPARENT VALUE OF MAGNITUDE
        if mag < m_min_current:
            continue
        nr_eq_current = nr_eq_current + 1  # NR OF EQ - s FROM m_min_current
        mean_mag += mag
        earthquakes.append({"magnitude": mag})
        if mag > m_max_obs:
            m_max_obs = mag
    # configuration['m_min'] = config_m_min
    if nr_eq_current > 0:
        mean_mag = mean_mag / nr_eq_current;
        b_aki_utsu = 1.0 / (np.log(10) * (mean_mag - m_min_current))
        print(f"ESTIMATED b-value (AKI-UTSU) = {b_aki_utsu}, beta = {b_aki_utsu * np.log(10)}")
    print(f"NUMBER of EQ-s (from M_min_total) = {nr_eq_total}")
    print(f"NUMBER of EQ-s (from M_min_current) = {nr_eq_current}")
    print(f"Maximum OBSERVED magnitude (APPARENT)  = {m_max_obs}")
    return {'begin': begin_time, 'end': end_time, 'time_span': time_span, 'm_min': m_min_current,
            'sd': sd_mag, 'name': name, 'earthquakes': earthquakes}, m_max_obs, m_min_current
    # END OF COMPLETE DATA COMPUTATIONS ====================================


def synthetic_catalogs_generation(configuration):
    """
    The procedure generates pseudo paleo, historical, and complete synthetic catalogues
    based on magnitude occurrence probability definition
    and catalogs configurations defined in the simulation part of the configuration.
    *Warning* The older catalogs are replaced by synthetic ones
    and new minimum and maximum observed magnitudes of all catalogues are set.

    :param configuration: The dictionary of all Ha3Py parameters.
        The function doesn't change the input configuration.
    :type configuration: dict
    :return: The new configuration with synthetic catalogues
    """
    # ----
    # Init
    # ----
    output_config = configuration.copy()
    simulation_config = configuration.get('simulation')
    if simulation_config is None:
        print('Missing the simulation definition in the configuration')
        return
    prehistoric_conf = simulation_config.get('paleo_catalog')
    historic_conf = simulation_config.get('historic_catalog')
    complete_conf = simulation_config.get('complete_catalogs')
    general_m_min = 100.0
    general_m_max_obs = -100.0
    general_sd_m_max_obs = 0.0
    # -------------
    # Paleo catalog
    # -------------
    if prehistoric_conf is None:
        output_config.pop('paleo_catalog', None)
    else:
        catalog, m_max_obs, m_min = prehistorical_simulation(prehistoric_conf, configuration)
        if m_max_obs > general_m_max_obs:
            general_m_max_obs = m_max_obs
            general_sd_m_max_obs = catalog['sd']
        if general_m_min > m_min:
            general_m_min = m_min
        if catalog['earthquakes']:
            output_config['paleo_catalog'] = catalog
        else:
            output_config.pop('paleo_catalog', None)
    # ----------------
    # Historic catalog
    # ----------------
    if historic_conf is None:
        output_config.pop('historic_catalog', None)
    else:
        catalog, m_max_obs, m_min = historical_simulation(historic_conf, configuration)
        if m_max_obs > general_m_max_obs:
            general_m_max_obs = m_max_obs
            general_sd_m_max_obs = catalog['sd']
        if general_m_min > m_min:
            general_m_min = m_min
        if catalog['earthquakes']:
            output_config['historic_catalog'] = catalog
        else:
            output_config.pop('historic_catalog', None)
    # -----------------
    # Complete catalogs
    # -----------------
    if complete_conf is None:
        output_config.pop('complete_catalogs', None)
    else:
        catalogs = []
        output_config.pop('complete_catalogs', None)
        for idx, catalog_conf in enumerate(complete_conf):
            catalog, m_max_obs, m_min = complete_simulation(catalog_conf, configuration, idx)
            if m_max_obs > general_m_max_obs:
                general_m_max_obs = m_max_obs
                general_sd_m_max_obs = catalog['sd']
            if general_m_min > m_min:
                general_m_min = m_min
            if catalog['earthquakes']:
                # test_lambda_beta(catalog)
                catalogs.append(catalog)
            if m_max_obs > general_m_max_obs:
                general_m_max_obs = m_max_obs
                general_sd_m_max_obs = catalog['sd']
            if general_m_min > m_min:
                general_m_min = m_min
        if catalogs:
            output_config['complete_catalogs'] = catalogs
        else:
            output_config.pop('complete_catalogs', None)
    # --------------
    # Final settings
    # --------------
    if general_m_max_obs > -99.0:
        output_config['m_max_obs'] = general_m_max_obs
        output_config['sd_m_max_obs'] = general_sd_m_max_obs
    else:
        output_config.pop('m_max_obs', None)
        output_config.pop('sd_m_max_obs', None)
    if general_m_min < 99.0:
        output_config['m_min'] = general_m_min
    else:
        output_config.pop('m_min', None)
    output_config.pop('m_max_current', None)
    # ---
    # End
    # ---
    return output_config


def main():
    configuration = load_configuration()
    output = synthetic_catalogs_generation(configuration)
    save_configuration(output)


if __name__ == "__main__":
    main()