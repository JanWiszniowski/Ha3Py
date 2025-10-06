from scipy.stats import truncnorm, uniform
from ha3py.utils import HaPyException
from ha3py.get_events_occurrence import get_events_occurrence
from ha3py.configuration import load_configuration, save_configuration
import numpy as np
import scipy.integrate as integrate


def time_cdf(t, m, occurrence_probability):
    r"""
    The cumulative distribution function of non occurrence event probability given magnitude.
    If there exists time_cdf method in the event_occurrence object,
    the result of time_cdf is returned, else

    .. math::
        F_T(t|m)=\frac{\int_{0}^{t} f_{M}(m,t)dt}{\int_{0}^{\infty } f_{M}(m,t)dt}

    is returned.

    :param t: The time
    :type t: float
    :param m: The magnitude
    :type m: float
    :param occurrence_probability: The object describing the event occurrence probability.
    :type occurrence_probability: OccurrenceBase
    :return: The CDF of the probability of the non occurrence of the event
        with the magnitude :math:`m` or greater in the time :math:`t`.
    :rtype: float
    """
    if callable(getattr(occurrence_probability, 'time_cdf', None)):
        return occurrence_probability.time_cdf(t, m)
    nominator = integrate.quad(lambda x: occurrence_probability.pdf(m, x), 0, t)
    denominator = integrate.quad(lambda x: occurrence_probability.pdf(m, x), 0, np.inf)
    return nominator[0] / denominator[0]


def time_rg(occurrence_probability, m=None):
    """
    Next event time random generator.
    It generates random period of non occurrence of events with magnitude grater than selected.
    If none magnitude is set, the minimum magnitude of the event occurrence object is taken.

    :param occurrence_probability: The object describing the event occurrence probability.
    :type occurrence_probability: OccurrenceBase
    :param m: The magnitude of non occurrence
    :type m: float
    :return: The random value of non occurrence period in years
    :rtype: float
    """
    if m is None:
        m = occurrence_probability.m_min
    urv = uniform.rvs()
    lower_limit = 0.0
    upper_limit = 1.0
    while time_cdf(upper_limit, m, occurrence_probability) < urv:
        upper_limit *= 2.0
    while upper_limit-lower_limit > 1e-8:
        new_point = (lower_limit+upper_limit) / 2.0
        value = time_cdf(new_point, m, occurrence_probability)
        if value > urv:
            upper_limit = new_point
        else:
            lower_limit = new_point
    return (lower_limit+upper_limit) / 2.0


class MagnitudeRandomise:
    """
    TODO The randomization of magnitudes does not work correctly.

    Initialisation parameters:

        *occurrence_probability:* (OccurrenceBase) The object describing the event occurrence probability.

        *magnitude_uncertainty:* (float) The standard deviation of simulated magnitude uncertainty

    """
    def __init__(self, occurrence_probability, magnitude_uncertainty):
        """
        Initialisation

        :param occurrence_probability: The object describing the event occurrence probability.
        :type occurrence_probability: OccurrenceBase
        :param magnitude_uncertainty: The standard deviation of simulated magnitude uncertainty
        :type magnitude_uncertainty: (float)
        """
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
        """
        Return random magnitude

        :return: Random magnitude
        :rtype: float
        """
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


def test_lambda_beta(catalog_configuration):
    r"""
    The function estimated and prints :math:`\beta` and :math:`\lambda` parameters
    of the catalogue assessed by the simplest method:

    .. math::
        \beta =\frac{1}{n}\sum_{i=1}^{n}m_i -m_{min},

    and

    .. math::
        \lambda = \frac{n}{T},

    where :math:`T` id the catalogue time span, :math:`n` is the number of events,
    and :math:`m_{min}` the minimum magnitude of completeness.

    :param catalog_configuration: The catalog definition, part of the dictionary of all Ha3Py parameters
    :type catalog_configuration: dict
    :return: None
    """
    m_min = catalog_configuration['m_min']
    magnitudes = [earthquake['magnitude'] for earthquake in catalog_configuration['earthquakes']]
    time_span = catalog_configuration['time_span']
    no_earthquakes = len(magnitudes)
    beta_for_catalog = 1.0 / (np.mean(magnitudes) - m_min)
    # magnitude_distribution = GutenbergRichter(configuration, beta=beta_for_catalog)
    lambda_for_catalog = float(no_earthquakes) / time_span
    print(f"Catalog {catalog_configuration['name']}: beta={beta_for_catalog}, lambda={lambda_for_catalog} for m_min={m_min}")


def extreme_simulation(catalog_configuration, occurrence_probability, name):
    r"""
    The function creates the synthetic extreme earthquakes.
    It simulates the historical or paleo catalogs, where data are incomplete.
    The time span of te catalogue is divided into random periods.
    In each period the random number of event with random magnitudes are generated
    and the maximum magnitude are added to the catalogue.
    The randomisation depends on the predefined occurrence probability.

    :param catalog_configuration: The catalogue definition, part of the dictionary of all Ha3Py parameters.
    :type catalog_configuration: dict
    :param occurrence_probability: The object describing the event occurrence probability.
    :type occurrence_probability: OccurrenceBase
    :param name: The name of the created synthetic catalogue.
    :type name: str
    :return: The synthetic catalogue
    :rtype: dict

    The simulation catalog configuration dictionary should contain following fields:

    **time_interval**: (float) The time interval,
    from which the maximum magnitude are put to the synthetic catalogue.

    **time_uncertainty**: (float of uniform probability) The uncertainty of periods time.

    **magnitude_uncertainty**: (float) The uncertainty for magnitude simulation.

    **sd**: (float) The standard deviation of catalog magnitudes.
    Required for compatibility, when magnitude_uncertainty is not defined.
    At last magnitude_uncertainty or sd must be defined.

    **m_min**: (str) The minimum completeness magnitude of the simulated catalog.
    Event with lower magnitude are removed.

    **begin, end, time_span**: The beginning, end, and time span of the catalog.
    At least two of them are required.

    """
    time_uncertainty = catalog_configuration.get('time_uncertainty', 10.0) / 2.0
    time_interval = catalog_configuration.get('time_interval', 20.0)
    begin_time, end_time, time_span = get_times(catalog_configuration)
    magnitude_uncertainty = catalog_configuration.get('magnitude_uncertainty')
    sd = catalog_configuration.get('sd', magnitude_uncertainty)
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
    first the non occurrence time of an event with magnitude greater the :math:`m_{min}` is generated,
    then the event time is set, next the magnitude is random generated for that event,
    and next events times and magnitudes are generated until reaching the end of the interval time.

    :param catalog_configuration: The catalogue definition, part of the dictionary of all Ha3Py parameters.
    :type catalog_configuration: dict
    :param occurrence_probability: The object describing the event occurrence probability.
    :type occurrence_probability: OccurrenceBase
    :param name: The name of the created synthetic catalogue.
    :type name: str
    :return: The synthetic catalogue
    :rtype: dict
    
    The simulation catalog configuration dictionary should contain following fields:

    **magnitude_uncertainty**: (float) The uncertainty for magnitude simulation.

    **sd**: (float) The standard deviation of catalog magnitudes.
    Required for compatibility, when magnitude_uncertainty is not defined.
    At last magnitude_uncertainty or sd must be defined.

    **m_min**: (str) The minimum completeness magnitude of the simulated catalog.
    Event with lower magnitude are removed.

    **begin, end, time_span**: The beginning, end, and time span of the catalog.
    At least two of them are required.

    """
    begin_time, end_time, time_span = get_times(catalog_configuration)
    magnitude_uncertainty = catalog_configuration.get('magnitude_uncertainty')
    sd = catalog_configuration.get('sd', magnitude_uncertainty)
    m_min = catalog_configuration['m_min']
    occurrence_probability.m_min = m_min
    earthquakes = []
    m_max_obs = -100.0
    dt = time_rg(occurrence_probability) * uniform.rvs()
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
        dt = time_rg(occurrence_probability)
        time += dt
    return {'begin': begin_time, 'end': end_time, 'time_span': time_span, 'm_min': m_min,
            'sd': sd, 'name': name, 'earthquakes': earthquakes}, m_max_obs, m_min


def full_simulation_without_date(catalog_configuration, occurrence_probability, name):
    r"""
    The full_simulation_none_date simulate the catalogue first by generation of random number events
    in the defined time interval based on the defined n-events occurrence probability
    and next for each event the magnitude is random generated based on the defined magnitude distribution.

    :param catalog_configuration: The catalogue definition, part of the dictionary of all Ha3Py parameters.
    :type catalog_configuration: dict
    :param occurrence_probability: The object describing the event occurrence probability.
    :type occurrence_probability: OccurrenceBase
    :param name: The name of the created synthetic catalogue.
    :type name: str
    :return: The synthetic catalogue
    :rtype: dict
    
    The simulation catalog configuration dictionary should contain following fields:

    **magnitude_uncertainty**: (float) The uncertainty for magnitude simulation.

    **sd**: (float) The standard deviation of catalog magnitudes.
    Required for compatibility, when magnitude_uncertainty is not defined.
    At last magnitude_uncertainty or sd must be defined.

    **m_min**: (str) The minimum completeness magnitude of the simulated catalog.
    Event with lower magnitude are removed.

    **begin, end, time_span**: The beginning, end, and time span of the catalog.
    At least two of them are required.

    """
    begin_time, end_time, time_span = get_times(catalog_configuration)
    magnitude_uncertainty = catalog_configuration.get('magnitude_uncertainty')
    sd = catalog_configuration.get('sd', magnitude_uncertainty)
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
    """
    The procedure remains the catalogue almost unchanged,
    but finds the minimum magnitude, if not defined, and maximum observed in the catalogue magnitude.
    If m_min is defined, removes event with magnitude smaller than m_min.
    The procedure doesn't estimate the catalogue completeness magnitude.

    :param catalog: The full catalogue description with earthquakes - not only simulation definition.
    :type catalog: dict
    :return: The copy of the input catalog with event below m_min removed,
        maximum observed in the catalogue magnitude, and minimum magnitude in the catalogue.
    :rtype: tuple
    """
    catalog = catalog.copy()
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
    occurrence_probability = get_events_occurrence(configuration)
    general_m_min = 100.0
    general_m_max_obs = -100.0
    general_sd_m_max_obs = 0.0
    # -------------
    # Paleo catalog
    # -------------
    if prehistoric_conf is None:
        output_config.pop('paleo_catalog', None)
    else:
        generator = prehistoric_conf.get('generator', 'full_simulation_incremental')
        name = prehistoric_conf.get('name', 'paleo_synthetic')
        if generator == 'full_simulation_incremental':
            catalog, m_max_obs, m_min = full_simulation_incremental(prehistoric_conf,
                                                                    occurrence_probability, name)
        elif generator == 'extreme_simulation':
            catalog, m_max_obs, m_min = extreme_simulation(prehistoric_conf,
                                                           occurrence_probability, name)
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
    # ----------------
    # Historic catalog
    # ----------------
    if historic_conf is None:
        output_config.pop('historic_catalog', None)
    else:
        generator = historic_conf.get('generator', 'extreme_simulation')
        name = historic_conf.get('name', 'historic_synthetic')
        if generator == 'full_simulation_incremental':
            catalog, m_max_obs, m_min = full_simulation_incremental(historic_conf,
                                                                    occurrence_probability, name)
        elif generator == 'extreme_simulation':
            catalog, m_max_obs, m_min = extreme_simulation(historic_conf,
                                                           occurrence_probability, name)
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
    # -----------------
    # Complete catalogs
    # -----------------
    if complete_conf is None:
        output_config.pop('complete_catalogs', None)
    else:
        catalogs = []
        old_catalogs = output_config.pop('complete_catalogs', None)
        for idx, catalog_conf in enumerate(complete_conf):
            generator = catalog_conf.get('generator', 'no_date_simulation')
            name = catalog_conf.get('name', f'complete#{idx}_synthetic')
            if generator == 'full_simulation_incremental':
                catalog, m_max_obs, m_min = full_simulation_incremental(catalog_conf,
                                                                        occurrence_probability, name)
            elif generator == 'no_date_simulation':
                catalog, m_max_obs, m_min = full_simulation_without_date(catalog_conf,
                                                                         occurrence_probability, name)
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