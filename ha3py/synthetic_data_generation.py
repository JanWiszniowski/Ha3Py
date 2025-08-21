from scipy.stats import truncnorm, uniform
from ha3py.utils import HaPyException
from ha3py.get_events_occurrence import get_events_occurrence
from ha3py.configuration import load_configuration, save_configuration

class MagnitudeRandomise:
    def __init__(self, occurrence_probability, magnitude_uncertainty):
        self.occurrence_probability = occurrence_probability
        self.magnitude_uncertainty = magnitude_uncertainty
        self.mu_3 = 3 * magnitude_uncertainty
        self.m_max = occurrence_probability.m_max
        self.m_max_bottom = self.m_max - self.mu_3

    def __call__(self):
        m = self.occurrence_probability.magnitude_distribution.rvs()
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
    time_uncertainty = catalog_configuration.get('time_uncertainty', 0.0) / 2.0
    time_interval = catalog_configuration.get('time_interval', 100.0)
    begin_time, end_time, time_span = get_times(catalog_configuration)
    magnitude_uncertainty = catalog_configuration.get('magnitude_uncertainty', 0.2)
    name = catalog_configuration.get('name', name)
    m_min = catalog_configuration['m_min']
    earthquakes = []
    m_max_obs = -100.0
    occurrence_probability.m_min = m_min - 3.0
    loop = True
    interval_end = begin_time
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
        event_time = truncnorm.rvs(interval_begin, interval_end, scale=interval/3)
        earthquakes.append({'magnitude': m_max_interval, 'date': event_time, 'time_span': interval,
                     'sd': magnitude_uncertainty})
        if m_max_obs < m_max_interval:
            m_max_obs = m_max_interval
    return {'begin': begin_time, 'end': end_time, 'time_span': time_span, 'm_min': m_min,
            'sd': magnitude_uncertainty, 'name': name, 'earthquakes': earthquakes}, m_max_obs, m_min


def full_simulation_incremental(catalog_configuration, occurrence_probability, name):
    r"""
    The incremental catalog generation: it starts with the beginning of the interval time,
    first the random non occurrence time, then the
    of event with magnitude greater the :math:`m_{min}` is generated, then the event time is set,
    next the magnitude is random generated for that event, and times and magnitudes are generated until
    reaching the end of the interval time.

    :param catalog_configuration: The catalog definition, part of the dictionary of all Ha3Py parameters
    :type catalog_configuration: dict
    :param occurrence_probability:
    :type occurrence_probability:
    :param name:
    :type name:
    :return:
    :rtype:
    """
    # time_uncertainty = catalog_configuration.get('time_uncertainty', 0.0)
    begin_time, end_time, time_span = get_times(catalog_configuration)
    magnitude_uncertainty = catalog_configuration.get('magnitude_uncertainty', 0.2)
    name = catalog_configuration.get('name', name)
    m_min = catalog_configuration['m_min']
    occurrence_probability.m_min = m_min - 3.0
    earthquakes = []
    m_max_obs = -100.0
    dt = occurrence_probability.rvs() * uniform.rvs()
    time = begin_time + dt
    magnitude_randomise = MagnitudeRandomise(occurrence_probability, magnitude_uncertainty)
    while time <= end_time:
        m = magnitude_randomise()
        if m > m_min:
            earthquakes.append({'magnitude': m, 'date': time, 'time_span': dt, 'sd': magnitude_uncertainty})
            if m > m_max_obs:
                m_max_obs = m
        dt = occurrence_probability.rvs()
        time += dt
    return {'begin': begin_time, 'end': end_time, 'time_span': time_span, 'm_min': m_min,
            'sd': magnitude_uncertainty, 'name': name, 'earthquakes': earthquakes}, m_max_obs, m_min


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
    magnitude_uncertainty = catalog_configuration.get('magnitude_uncertainty', 0.2)
    name = catalog_configuration.get('name', name)
    m_min = catalog_configuration['m_min']
    occurrence_probability.m_min = m_min - 3.0
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
            'sd': magnitude_uncertainty, 'name': name, 'earthquakes': earthquakes}, m_max_obs, m_min

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
    general_m_min = 100.0
    general_m_max_obs = -100.0
    general_sd_m_max_obs = 0.0
    # Paleo catalog
    if prehistoric_conf is None:
        output_config.pop('paleo_catalog', None)
    else:
        generator = prehistoric_conf.get('generator', 'full_simulation_incremental')
        if generator == 'full_simulation_incremental':
            catalog, m_max_obs, m_min = full_simulation_incremental(prehistoric_conf, event_occurrence, 'paleo')
        elif generator == 'extreme_simulation':
            catalog, m_max_obs, m_min = extreme_simulation(prehistoric_conf, event_occurrence, 'paleo')
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
        if generator == 'full_simulation_incremental':
            catalog, m_max_obs, m_min = full_simulation_incremental(historic_conf, event_occurrence, 'historic')
        elif generator == 'extreme_simulation':
            catalog, m_max_obs, m_min = extreme_simulation(historic_conf, event_occurrence, 'historic')
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
        old_catalogs = output_config.pop('complete_catalogs')
        for idx, catalog_conf in enumerate(complete_conf):
            generator = catalog_conf.get('generator', 'full_simulation_incremental')
            if generator == 'full_simulation_incremental':
                catalog, m_max_obs, m_min = full_simulation_incremental(catalog_conf, event_occurrence, f'complete#{idx}')
            elif generator == 'no_date_simulation':
                catalog, m_max_obs, m_min = full_simulation_without_date(catalog_conf, event_occurrence, f'complete#{idx}')
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
                output_config['historic_catalog'] = catalog
            else:
                output_config.pop('historic_catalog', None)
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