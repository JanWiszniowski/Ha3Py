from scipy.stats import truncnorm, uniform
from ha3py.utils import HaPyException
from ha3py.get_events_occurrence import get_events_occurrence
from ha3py.configuration import load_configuration, save_configuration
import numpy as np
import scipy.integrate as integrate
from scipy.optimize import fsolve


def time_cdf(t, m, event_occurrence):
    if callable(getattr(event_occurrence, 'time_cdf', None)):
        tcdf = event_occurrence.time_cdf(t, m)
        return tcdf
    nominator = integrate.quad(lambda x: event_occurrence.pdf(m, x), 0, t)
    denominator = integrate.quad(lambda x: event_occurrence.pdf(m, x), 0, np.inf)
    return nominator[0] / denominator[0]

# def _fun(t, event_occurrence, m, urv ):
#     if t < 0:
#         f = 1.0e6 * (np.exp(-t))
#     else:
#     tcdf = time_cdf(t, m, event_occurrence)
#     f = tcdf - urv
#     f += 1.0e6 * t ** 8 if t < 0 else 0
#     print(f"urv {t}, {tcdf}, {urv}, {f}")
#     return f
#
# def time_rng(event_occurrence, m=None):
#     if m is None:
#         m = event_occurrence.m_min
#     urv = uniform.rvs()
#     print(f"urv {urv}")
#     time_rv = fsolve(_fun, np.array(1.0), args=(event_occurrence, m, urv))
#     return time_rv[0]

def time_rng(event_occurrence, m=None):
    if m is None:
        m = event_occurrence.m_min
    urv = uniform.rvs()
    # print(f"urv={urv}")
    lower_limit = 0.0
    upper_limit = 1.0
    while time_cdf(upper_limit, m, event_occurrence) < urv:
        upper_limit *= 2.0
    while upper_limit-lower_limit > 1e-8:
        new_point = (lower_limit+upper_limit) / 2.0
        value = time_cdf(new_point, m, event_occurrence)
        # print(f"==>{lower_limit} {upper_limit} {new_point} {value-urv}")
        if value > urv:
            upper_limit = new_point
        else:
            lower_limit = new_point
    return (lower_limit+upper_limit) / 2.0


class MagnitudeRandomise:
    def __init__(self, occurrence_probability, magnitude_uncertainty):
        self.occurrence_probability = occurrence_probability
        if magnitude_uncertainty is None:
            self.magnitude_uncertainty = None
            self.mu_3 = None
            self.m_max = None
            self.m_max_bottom = None
        else:
            self.occurrence_probability.m_min = occurrence_probability.m_min - 4.0 * magnitude_uncertainty
            self.magnitude_uncertainty = magnitude_uncertainty
            self.mu_3 = 3.0 * magnitude_uncertainty
            self.m_max = occurrence_probability.m_max
            self.m_max_bottom = self.m_max - self.mu_3

    def __call__(self):
        m = self.occurrence_probability.magnitude_distribution.rvs()
        if self.magnitude_uncertainty is None:
            return m
        if m < self.m_max_bottom:
            m += truncnorm.rvs(-self.mu_3, self.mu_3, scale=self.magnitude_uncertainty)
        else:
            m += truncnorm.rvs(-self.mu_3, self.m_max - m, scale=self.magnitude_uncertainty)
        return m

def get_times(catalog_configuration):
    """

    :param catalog_configuration: The catalog definition, part of the dictionary of all Ha3Py parameters
    :type catalog_configuration:
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

def test_lambda_beta(catalog):
    m_min = catalog['m_min']
    magnitudes = [earthquake['magnitude'] for earthquake in catalog['earthquakes']]
    time_span = catalog['time_span']
    no_earthquakes = len(magnitudes)
    beta_for_catalog = 1.0 / (np.mean(magnitudes) - m_min)
    # magnitude_distribution = GutenbergRichter(configuration, beta=beta_for_catalog)
    lambda_for_catalog = float(no_earthquakes) / time_span
    print(f"Catalog {catalog['name']}: beta={beta_for_catalog}, lambda={lambda_for_catalog} for m_min={m_min}")


def extreme_simulation(catalog_configuration, occurrence_probability, name):
    r"""
    Extreme simulation simulate extreme earthquakes
    :param catalog_configuration: The catalog definition, part of the dictionary of all Ha3Py parameters.
    :type catalog_configuration: dict
    :param occurrence_probability:
    :type occurrence_probability:
    :param name: The name of the catalogue.
    :type name: (str)
    :return: The catalogue
    :rtype: dict

    The catalog configuration dictionary should contain following fields:
    :time_uncertainty:
    :time_interval:
    :magnitude_uncertainty: The uncertainty for magnitude simulation (required)
    :name: The name of the catalog
    :m_min: The minimum completeness magnitude of the simulated catalog. Event with lower magnitude are removed
    :begin, end, time_span: The beginning, end, and time span of the catalog. At least two of them are required.
    """
    time_uncertainty = catalog_configuration.get('time_uncertainty', 10.0) / 2.0
    time_interval = catalog_configuration.get('time_interval', 10.0)
    begin_time, end_time, time_span = get_times(catalog_configuration)
    magnitude_uncertainty = catalog_configuration.get('magnitude_uncertainty', 0.2)
    sd = catalog_configuration.get('sd', magnitude_uncertainty)
    name = catalog_configuration.get('name', name)
    m_min = catalog_configuration['m_min']
    earthquakes = []
    m_max_obs = -100.0
    occurrence_probability.m_min = m_min
    loop = True
    interval_end = begin_time
    event_time = begin_time
    magnitude_randomise = MagnitudeRandomise(occurrence_probability, magnitude_uncertainty)
    while loop:
        interval_begin = interval_end
        interval_end = interval_begin + time_interval + uniform.rvs(-time_uncertainty, time_uncertainty)
        if interval_end > end_time:
            interval_end = end_time
            loop = False
        interval = interval_end - interval_begin
        n = occurrence_probability.d_rvs(interval)
        m_max_interval = -100.0
        if n == 0:
            continue
        for idx in range(n):
            m = magnitude_randomise()
            if m_max_interval < m:
                m_max_interval = m
        if m_max_interval < m_min:
            continue
        last_event_time = event_time
        event_time = (interval_begin + interval_end) / 2.0 + truncnorm.rvs(-interval/2, interval/2,
                                                                          scale=interval/6)
        earthquakes.append({'magnitude': m_max_interval, 'date': event_time,
                            'time_span': event_time-last_event_time, 'sd': sd})
        if m_max_obs < m_max_interval:
            m_max_obs = m_max_interval
    return {'begin': begin_time, 'end': end_time, 'time_span': time_span, 'm_min': m_min,
            'sd': sd, 'name': name, 'earthquakes': earthquakes}, m_max_obs, m_min


def full_simulation_incremental(catalog_configuration, occurrence_probability, name):
    r"""
    The incremental catalog generation: it starts with the beginning of the interval time,
    first the non occurrence time  of an event with magnitude greater the :math:`m_{min}` is generated,
    then the event time is set, next the magnitude is random generated for that event,
    and next events times and magnitudes are generated until reaching the end of the interval time.

    :param catalog_configuration: The catalog definition, part of the dictionary of all Ha3Py parameters
    :type catalog_configuration: dict
    :param occurrence_probability:
    :type occurrence_probability:
    :param name: name of the simulated catalogue
    :type name: str
    :return:
    :rtype:
    """
    # time_uncertainty = catalog_configuration.get('time_uncertainty', 0.0)
    begin_time, end_time, time_span = get_times(catalog_configuration)
    magnitude_uncertainty = catalog_configuration.get('magnitude_uncertainty', 0.2)
    sd = catalog_configuration.get('sd', magnitude_uncertainty)
    # time_step = catalog_configuration.get('time_step', 1.0)
    name = catalog_configuration.get('name', name)
    m_min = catalog_configuration['m_min']
    occurrence_probability.m_min = m_min
    earthquakes = []
    m_max_obs = -100.0
    dt = time_rng(occurrence_probability) * uniform.rvs()
    # dt = occurrence_probability.t_rng() * uniform.rvs()
    # m = occurrence_probability.rvs(dt) * (1.0 + magnitude_uncertainty * (uniform.rvs() -0.5))
    time = begin_time + dt
    magnitude_randomise = MagnitudeRandomise(occurrence_probability, magnitude_uncertainty)
    while time <= end_time:
        m = magnitude_randomise()
        if m > m_min:
            earthquakes.append({'magnitude': m, 'date': time, 'time_span': dt, 'sd': sd})
            if m > m_max_obs:
                m_max_obs = m
        # dt = occurrence_probability.t_rng()
        # m = occurrence_probability.rvs(dt) * (1.0 + magnitude_uncertainty * (uniform.rvs() - 0.5))
        dt = time_rng(occurrence_probability)
        time += dt
    return {'begin': begin_time, 'end': end_time, 'time_span': time_span, 'm_min': m_min,
            'sd': sd, 'name': name, 'earthquakes': earthquakes}, m_max_obs, m_min


def full_simulation_without_date(catalog_configuration, occurrence_probability, name):
    r"""
    The full_simulation_none_date simulate the catalogue first by generation of random number events
    in the defined time interval based on the defined n-events occurrence probability
    and next for each event the magnitude is random generated based on the defined magnitude distribution.

    :param catalog_configuration: The catalog definition, part of the dictionary of all Ha3Py parameters
    :type catalog_configuration: dict
    :param occurrence_probability:
    :type occurrence_probability:
    :param name:
    :type name:
    :return:
    :rtype:
    """
    begin_time, end_time, time_span = get_times(catalog_configuration)
    magnitude_uncertainty = catalog_configuration.get('magnitude_uncertainty')
    sd = catalog_configuration.get('sd', magnitude_uncertainty)
    name = catalog_configuration.get('name', name)
    m_min = catalog_configuration['m_min']
    occurrence_probability.m_min = m_min
    n = occurrence_probability.d_rvs(time_span)
    if n == 0:
        return None
    earthquakes = []
    m_max_obs = -100.0
    magnitude_randomise = MagnitudeRandomise(occurrence_probability, magnitude_uncertainty)
    for idx in range(n):
        m = magnitude_randomise()
        if m > m_min:
            earthquakes.append({'magnitude': m})
            if m > m_max_obs:
                m_max_obs = m
    return {'begin': begin_time, 'end': end_time, 'time_span': time_span, 'm_min': m_min,
            'sd': sd, 'name': name, 'earthquakes': earthquakes}, m_max_obs, m_min

def remain_unchanged(catalog):
    earthquakes = catalog.get('earthquakes')
    if earthquakes is None:
        print(f"Incorrect catalog {catalog.get('name', 'without name')} won't remain")
        return None, -100.0, 100.0
    m_min = catalog.get('m_min')
    m_max_obs = -100.0
    if m_min in None:
        m_min = 100.0
        for earthquake in earthquakes:
            if m_min > earthquake.magnitude:
                m_min = earthquake.magnitude
        catalog['m_min'] = m_min
    else:
        earthquakes = [earthquake for earthquake in earthquakes  if earthquake.magnitude >= m_min]
        catalog['earthquakes'] = earthquakes
    for earthquake in earthquakes:
        if m_max_obs < earthquake.magnitude:
            m_max_obs = earthquake.magnitude
    return catalog, m_max_obs, m_min


def data_generation(configuration):
    """
    The procedure generates the synthetic catalogue based on magnitude occurrence probability definition
    and catalogs configurations defined in the simulation part of configuration.
    The older catalogs are replaced by synthetic ones.

    :param configuration: The dictionary of all Ha3Py parameters
    :type configuration: dict
    :return:
    """
    # Init
    output_config = configuration.copy()
    simulation_config = configuration.get('simulation')
    if simulation_config is None:
        print('Missing the simulation definition in the configuration')
        return
    prehistoric_conf = simulation_config.get('paleo_catalog')
    historic_conf = simulation_config.get('historic_catalog')
    complete_conf = simulation_config.get('complete_catalogs')
    event_occurrence = get_events_occurrence(configuration)
    # m_max_sim = 0
    # idx = 0
    # event_occurrence.m_min = event_occurrence.m_min - 1.0
    # while True:
    #     m = event_occurrence.magnitude_distribution.rvs()
    #
    #     if m_max_sim < m:
    #         m_max_sim = m
    #     print(f"{idx}: {m_max_sim} > {m}")
    #     idx +=1
    general_m_min = 100.0
    general_m_max_obs = -100.0
    general_sd_m_max_obs = 0.0
    # Paleo catalog
    if prehistoric_conf is None:
        output_config.pop('paleo_catalog', None)
    else:
        generator = prehistoric_conf.get('generator', 'full_simulation_incremental')
        name = prehistoric_conf.get('name', 'paleo_synthetic')
        if generator == 'full_simulation_incremental':
            catalog, m_max_obs, m_min = full_simulation_incremental(prehistoric_conf, event_occurrence, name)
        elif generator == 'extreme_simulation':
            catalog, m_max_obs, m_min = extreme_simulation(prehistoric_conf, event_occurrence, name)
        elif generator == 'remain unchanged':
            catalog, m_max_obs, m_min = remain_unchanged(configuration.get('paleo_catalog'))
        else:
            raise HaPyException(f'Generator {generator} not available for the paleo catalog')
        if m_max_obs > general_m_max_obs:
            general_m_max_obs = m_max_obs
            general_sd_m_max_obs = catalog['sd']
        if general_m_min > m_min:
            general_m_min = m_min
        if catalog['earthquakes']:
            output_config['paleo_catalog'] = catalog
        else:
            output_config.pop('paleo_catalog', None)
    # Historic catalog
    if historic_conf is None:
        output_config.pop('historic_catalog', None)
    else:
        generator = historic_conf.get('generator', 'extreme_simulation')
        name = historic_conf.get('name', 'historic_synthetic')
        if generator == 'full_simulation_incremental':
            catalog, m_max_obs, m_min = full_simulation_incremental(historic_conf, event_occurrence, name)
        elif generator == 'extreme_simulation':
            catalog, m_max_obs, m_min = extreme_simulation(historic_conf, event_occurrence, name)
        elif generator == 'remain unchanged':
            catalog, m_max_obs, m_min = remain_unchanged(configuration.get('historic_catalog'))
        else:
            raise HaPyException(f'Generator {generator} not available for the historic catalog')
        if m_max_obs > general_m_max_obs:
            general_m_max_obs = m_max_obs
            general_sd_m_max_obs = catalog['sd']
        if general_m_min > m_min:
            general_m_min = m_min
        if catalog['earthquakes']:
            output_config['historic_catalog'] = catalog
        else:
            output_config.pop('historic_catalog', None)
    # Complete catalogs
    if complete_conf is None:
        output_config.pop('complete_catalogs', None)
    else:
        catalogs = []
        old_catalogs = output_config.pop('complete_catalogs', None)
        for idx, catalog_conf in enumerate(complete_conf):
            generator = catalog_conf.get('generator', 'no_date_simulation')
            name = catalog_conf.get('name', f'complete#{idx}_synthetic')
            if generator == 'full_simulation_incremental':
                catalog, m_max_obs, m_min = full_simulation_incremental(catalog_conf, event_occurrence, name)
            elif generator == 'no_date_simulation':
                catalog, m_max_obs, m_min = full_simulation_without_date(catalog_conf, event_occurrence, name)
            elif generator == 'remain unchanged':
                if old_catalogs is not None and len(old_catalogs) > idx:
                    catalog, m_max_obs, m_min = remain_unchanged(old_catalogs[idx])
                else:
                    catalog, m_max_obs, m_min = (None, -100.0, 100.0)
                    print(f"Complete catalog {idx} don't exist. It won't remain.")
            else:
                raise HaPyException(f'Generator {generator} not available for the complete#{idx} catalog')
            if m_max_obs > general_m_max_obs:
                general_m_max_obs = m_max_obs
                general_sd_m_max_obs = catalog['sd']
            if general_m_min > m_min:
                general_m_min = m_min
            if catalog['earthquakes']:
                test_lambda_beta(catalog)
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
    # Final settings
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
    # End
    return output_config


def main():
    configuration = load_configuration()
    output = data_generation(configuration)
    save_configuration(output)


if __name__ == "__main__":
    main()