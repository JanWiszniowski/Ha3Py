r"""
Import catalogues to the Ha3Py configuration
--------------------------------------------

..
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
import numpy as np
import scipy as sp
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.event.event import Event
from obspy.core.event.origin import Origin
from obspy.core.event.magnitude import Magnitude
from obspy.core.event import read_events
from obspy.core.event.base import Comment
from ha3py.configuration import date_years, load_configuration, save_configuration


class EPISODESCatalog:
    """
    The class EPISODESCatalog allows read the EPISODES catalogue.
    It keeps the catalogue events in the ObsPy format.
    Therefore, it is possible to operate on the EPISODESCatalog events similar to the ObsPy catalogue.

    """

    def __init__(self):
        self.events = []


def _get_vals(cat, field):
    fields = cat['field'].flatten()
    idx = np.where(fields == [field])
    if not np.size(idx):
        return []
    vals_tab = cat['val'].flatten()
    vals = vals_tab[idx]
    return vals[0]


def read_episodes(file_name):
    """
    It reads the catalogue Matlab file and create the EPISODES catalog.

    :param file_name:
    :type file_name: str
    :return: The catalog
    :rtype: EPISODESCatalog
    """
    catalog = EPISODESCatalog()
    contents = sp.io.loadmat(file_name)
    for key in contents:
        if not key[0:2] == '__':
            break
    else:
        raise 'Missing variable in the EPOS-AH file'
    cat = contents[key]
    indexes = _get_vals(cat, 'ID')
    times = _get_vals(cat, 'Time')
    latitudes = _get_vals(cat, 'Lat')
    longitudes = _get_vals(cat, 'Long')
    depths = _get_vals(cat, 'Depth')
    mls = _get_vals(cat, 'ML')
    mws = _get_vals(cat, 'Mw')
    energy = _get_vals(cat, 'E')
    events_size = indexes.size
    for idx in range(events_size):
        utc_time = UTCDateTime((times[idx][0] - 719529.0) * 86400.0)
        origin = Origin(time=utc_time, longitude=longitudes[idx][0], latitude=latitudes[idx][0], depth=depths[idx][0])
        magnitudes = []
        if np.any(mls) and np.size(mls[idx]) and not np.isnan(mls[idx]):
            magnitudes.append(Magnitude(mag=mls[idx], magnitude_type='ML'))
        if np.any(mws) and np.size(mws[idx]) and not np.isnan(mws[idx]):
            magnitudes.append(Magnitude(mag=mws[idx], magnitude_type='Mw'))
        if np.any(energy) and np.size(energy[idx]) and not np.isnan(energy[idx]):
            magnitudes.append(Magnitude(mag=energy[idx], magnitude_type='Energy'))
        event = Event(resource_id=indexes[idx][0][0], origins=[origin], magnitudes=magnitudes)
        catalog.events.append(event)
    return catalog


def get_magnitude(event, magnitude_type=None):
    """
    Function get_magnitude extracts the magnitude of the event.
    If you want to extract a specific magnitude, you can define it as magnitude_type,
    e.g. *get_magnitude(event, magnitude_type='Mw')*, otherwise, any magnitude will be extracted.
    If the preferred_magnitude_id of the event is set, it returns the preferred origin.
    Otherwise, it returns the first magnitude from the list.
    The function is intended to extract the magnitude unconditionally and non-interactively.
    Therefore, if preferred_magnitude_id is not set and there are multiple magnitudes,
    the returned origin may be random.

    If event magnitude does not exist, but station_name magnitudes exist, the new magnitude is computed
    as the mean value of station_name magnitudes.

    :param event: The event object
    :type event: ObsPy Event
    :param magnitude_type:  (optional)
        Describes the type of magnitude. This is a free text. Proposed values are:
            * unspecified magnitude ('M') - function search for exactly unspecified magnitude,
            * local magnitude ('ML'),
            * moment magnitude ('Mw'),
            * energy ('Energy'),
            * etc.
    :type magnitude_type: str

    :return: The magnitude object or None if the function cannot find or create the magnitude.
        If only station_name magnitudes exist, the new ObsPy Magnitude object is created,
        but it is not appended to the event
    :rtype: ObsPy Magnitude

    """
    if event.preferred_magnitude_id is not None:
        magnitude_object = event.preferred_magnitude_id.get_referred_object()
        if not magnitude_type or magnitude_object.magnitude_type == magnitude_type:
            return magnitude_object
    if not event.magnitudes:
        if event.station_magnitudes:
            no_magnitudes = 0
            magnitude = 0.0
            if magnitude_type:
                for m in event.station_magnitudes:
                    if m.magnitude_type == magnitude_type:
                        magnitude += m.mag
                        no_magnitudes += 1
            else:
                for m in event.station_magnitudes:
                    magnitude += m.mag
                    no_magnitudes += 1
            if no_magnitudes:
                magnitude /= no_magnitudes
                magnitude_object = Magnitude(mag=magnitude, magnitude_type=magnitude_type)
                magnitude_object.comment.append(Comment(text=f"Mean of {no_magnitudes} station_name magnitudes"))
    else:
        if magnitude_type:
            for m in event.magnitudes:
                if m.magnitude_type == magnitude_type:
                    return m
        else:
            return event.magnitudes[0]
    return None


def remove_catalogs(configuration):
    """
    The procedure lists complete catalogues, and allows an operator to select, which one to delete.

    :param configuration: The configuration of Ha3Py including catalogues.
    :type configuration: dict
    """
    catalogs = configuration.get('complete_catalogs')
    if catalogs is None:
        catalogs = []
        configuration['complete_catalogs'] = catalogs
    elif catalogs:
        loop = 'loop'
        while not loop and catalogs:
            print(f'There exists {len(catalogs)} complete catalogs:')
            for idx, catalog in enumerate(catalogs):
                print(f" {idx}: {catalog['name']}")
            loop = input("Write number of name of the catalog to remove it or ENTER to skip > ")
            if loop:
                if loop.isdigit():
                    del catalogs[int(loop)]
                else:
                    for idx, catalog in enumerate(catalogs):
                        if loop == catalog['name']:
                            del catalogs[idx]


def get_origin(event):
    """
    Function get_origin extracts the origin from the event.
    If the preferred_origin_id of the event is set, it returns the preferred origin.
    Otherwise, it returns the first origin from the list.
    The function is intended to extract the event origin unconditionally and non-interactively.
    Therefore, if preferred_origin_id is not set and there are multiple origins, the returned origin may be random

    :param event: The event object
    :type event: ObsPy Event

    :return: The origin object or None if no origin is defined for the event.
    :rtype: ObsPy Origin

    """
    if not event.origins:
        return None
    if event.preferred_origin_id is not None:
        return event.preferred_origin_id.get_referred_object()
    return event.origins[0]


def obspy_to_ha3py(obspy_catalog, magnitude_type='Mw'):
    """
    It converts the ObsPy or EPISODESCatalog to the dictionary object accepted by the Ha3Py configuration.

    :param obspy_catalog: The input catalog
    :type obspy_catalog: ObsPy Catalog or EPISODESCatalog
    :param magnitude_type: The name of magnitude (The default is recommended magnitude Mw)
    :type magnitude_type: str

    :return: The complete catalogue
    :rtype: dict
    """
    default_sd = float(input('Default standard deviation of magnitude > '))
    proposed_m_min = float(input('Minimum magnitude > '))
    catalog_m_min = 100.0
    catalog = dict()
    earthquakes = []
    sum_sd = 0.0
    begin_time = UTCDateTime(1000000, 1, 1, 0, 0)
    end_time = UTCDateTime(-1000000, 1, 1, 0, 0)
    for event in obspy_catalog.events:
        magnitude = get_magnitude(event, magnitude_type=magnitude_type)
        if magnitude is None:
            continue
        if magnitude.mag < proposed_m_min:
            continue
        origin = get_origin(event)
        if origin is None:
            continue
        earthquake = {'magnitude': magnitude.mag}
        if catalog_m_min > magnitude.mag:
            catalog_m_min = magnitude.mag
        if magnitude.mag_errors and magnitude.mag_errors.uncertainty:
            earthquake['sd'] = magnitude.mag_errors.uncertainty
        else:
            earthquake['sd'] = default_sd
        sum_sd += earthquake['sd']
        time = origin.time
        if begin_time > time:
            begin_time = time
        if end_time < time:
            end_time = time
        earthquake['date'] = date_years(time.year, time.month, time.day)
        earthquakes.append(earthquake)
    catalog['begin'] = date_years(begin_time.year, begin_time.month, begin_time.day)
    catalog['end'] = date_years(end_time.year, end_time.month, end_time.day + 1)
    catalog['time_span'] = catalog['end'] - catalog['begin']
    catalog['m_min'] = proposed_m_min
    catalog['sd'] = sum_sd / len(earthquakes)
    catalog['earthquakes'] = earthquakes
    catalog['name'] = obspy_catalog.resource_id.id
    return catalog


def read_file(file_name, catalog_format, configuration):
    """
    Procedure read the catalogue file.

    :param file_name: The path to the reading catalogue file
    :type file_name: str
    :param catalog_format: The catalogue format name. Allowed name is: EPISODES
    :type catalog_format: str
    :param configuration: The dictionary of all Ha3Py parameters
    :type configuration: dict
    :return:
    """
    if catalog_format == 'EPISODES':
        catalog = read_episodes(file_name)
    else:
        catalog = read_events(file_name, format=catalog_format)
    remove_catalogs(configuration)
    obspy_to_ha3py(catalog, configuration.get('Import magnitude type'))


def main():
    if len(sys.argv) < 4:
        print('Call import_to_hapy <configuration_file.json> <imported_catalogue> <format_name>')
        exit(0)
    params = load_configuration()
    read_file(sys.argv[2], sys.argv[2], params)
    save_configuration(params)


if __name__ == "__main__":
    main()
