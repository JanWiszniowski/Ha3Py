r"""
..
    Main program for assessment of seismic hazard parameters
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
from ha3py.compute import compute
from ha3py.print_results import print_results, print_hazard, print_percentage_share, compute_hazard
from ha3py.configuration import load_configuration, save_configuration, define_configuration
from ha3py.utils import HaPyException
from ha3py.plot_results import plot_results


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
    config_modified = define_configuration(configuration)
    if config_modified:
        choice_file = input('File to save the only configuration [enter if not save] > ')
        if choice_file:
            if choice_file.find(".") == -1:
                choice_file += '.json'
            save_configuration(configuration, file_name=choice_file)
    compute(configuration)
    save_configuration(configuration)
    print_percentage_share(configuration)
    print_results(configuration)
    configuration['hazard'] = compute_hazard(configuration)
    print_hazard(configuration)
    if 'output_text_file' in configuration:
        stdout_no = sys.stdout
        sys.stdout = open(configuration['output_text_file'], 'w')
        print_results(configuration)
        print_hazard(configuration)
        sys.stdout.close()
        sys.stdout = stdout_no
    plot_results(configuration)


if __name__ == "__main__":
    main()