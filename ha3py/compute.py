"""
Main module for probability of exceeding the magnitude estimation
-----------------------------------------------------------------

..
    :copyright:
        Jan Wiszniowski <jwisz@igf.edu.pl>,
        Andrzej Kijko <andrzej.kijko@up.ac.za>
    :license:
        GNU Lesser General Public License, Version 3
        (https://www.gnu.org/copyleft/lesser.html)
    :version 0.0.1:
        2025-02-01

"""

# import warnings
import datetime
import numpy as np
from scipy import optimize
from ha3py.ln_likelihood import ln_likelihood
from ha3py.constant_values import METHODS, LN_10_
from ha3py.compute_m_max import m_max_estimation
from ha3py.uncertainty import compute_uncertainty
from ha3py.get_events_occurrence import get_events_occurrence
from ha3py.corrections import lambda_correction
from ha3py.configuration import load_configuration, save_configuration, set_modify, input_configuration_name
from ha3py.print_results import print_percentage_share
from ha3py.utils import print_separation_double_line, print_separation_single_line


def opt_function(x_v, configuration, m_max=None):
    return ln_likelihood(np.exp(x_v), configuration, m_max=m_max)


def compute_occurrence(configuration):
    r"""
    The function estimates the events occurrence probability coefficients excluding m_max
    (e.g. :math:`\lambda` and :math:`\beta`)
    by finding the logarithm likelihood extremal value.

    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict
    :return:
        names of estimated coefficients,

        values of estimated coefficients
    :rtype: (tuple)

    """
    event_occurrence = get_events_occurrence(configuration)
    # x0_v = np.array(event_occurrence.coefficients[:-1])  # Array of beta and lambda for likelihood optimization
    log_x0_v = np.log(np.array(event_occurrence.coefficients[:-1]))  # Array of beta and lambda for likelihood optimization
    # results = optimize.minimize(ln_likelihood, x0_v, args=configuration,
    #                             method=configuration.get('likelihood_optimization_method'),
    #                             tol=configuration.get('likelihood_tolerance', 0.001))
    results = optimize.minimize(opt_function, log_x0_v, args=configuration,
                                method=configuration.get('likelihood_optimization_method'),
                                tol=configuration.get('likelihood_tolerance', 0.001))
    print('Likelihood estimation result: {}'.format(results.message))
    names = event_occurrence.coefficient_names[:-1]
    for idx, name in enumerate(names):
        configuration[name] = np.exp(results.x[idx])
        # configuration[name] = results.x[idx]
    # configuration['n'] = event_occurrence.d_expected()
    configuration['coefficient_names'] = event_occurrence.coefficient_names

    return names, np.exp(results.x)


def compute(configuration):
    r"""
    The procedure estimates the earthquake occurrence probability including magnitude distribution coefficients
    (e.g. estimates :math:`\lambda` and :math:`\beta`)
    and asses maximum event magnitude :math:`m_{max}` in the investigated area.
    The assessment method depends on the params.
    Results are saved to the output file. The existing in the params coefficients are replaced.

    The :math:`m_{max}` is assessed with the cooperation with the user.
    The procedure estimates and  proposes the :math:`m_{max}` value, which operator must confirm.

    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict
    :return: None

    """
    configuration.pop('mag_max', None)
    if 'm_max_current' not in configuration:
        configuration['m_max_current'] = configuration['m_max_obs'] + 0.5
    nr_runs = 0
    id_do_again = True
    suggested_m_max, sd_m_max = 9.0, 0.0

    while id_do_again:  # Loop of mag_max assessment
        nr_runs += 1
        names, values = compute_occurrence(configuration)
        lambda_correction(configuration)
        suggested_m_max, sd_m_max = m_max_estimation(configuration)
        configuration.pop('lambda_ref', None)
        print_separation_double_line()
        print(f"Run round number #{nr_runs}")
        print_separation_single_line()
        for name in names:
            if name == 'beta':
                print(f"{name:<10s} ={configuration[name]:7.3f} (b = {configuration[name]/LN_10_:3.1f})")
            elif name == 'lambda':
                print(f"{name:<10s} ={configuration[name]:7.3f} (for m_min = {configuration['m_min']:4.2f})")
            else:
                print(f"{name:<10s} ={configuration[name]:7.3f}")
        print(f"m_max (current) = {configuration['m_max_current']}")
        if suggested_m_max is None or suggested_m_max > 9.9:
            print(' I am SORRY. I DO NOT HAVE SUGGESTION REGARDING M_MAX')
        else:
            configuration['m_max_suggested'] = suggested_m_max
            print(f"Suggested value of m_max = {suggested_m_max:5.3f} (sd_m_max = {sd_m_max:5.3f})")
        print(f"                           for m_max_obs = {configuration['m_max_obs']:4.2f}")
        if configuration.get('bayesian_m_max_assessment'):
            print(f"                           for prior_m_max = {configuration['prior_m_max']:4.2f}")
        if configuration.get('enter_m_max', True):
            print_separation_double_line()
            if abs(suggested_m_max - configuration['m_max_current']) >= 0.01:
                print('To obtain the optimal solution for m_max, please re-run the process')
                print('according to the suggested value until the SUGGESTED and SOLUTION')
                print('values are the same.')
            else:
                print('')
                print('You have reach an optimal solution !!!')
            enter_correct = False
            print('')
            while not enter_correct:
                new_m_max_str = input('NEW value of m_max (NOT LESS than {:4.2f}) (or enter to accept current {:4.2f} and finish) >'.
                                      format(configuration['m_max_obs'], configuration['m_max_current']))
                if new_m_max_str:
                    try:
                        new_m_max = float(new_m_max_str)
                    except ValueError:
                        new_m_max = -1.0
                    if configuration['m_max_obs'] < new_m_max < 10.0:
                        configuration['m_max_current'] = new_m_max
                        enter_correct = True
                    else:
                        print('WRONG INPUT! "{}" is not a correct mag_max'.format(new_m_max_str))
                else:
                    id_do_again = False
                    enter_correct = True
        else:
            if suggested_m_max is not None and suggested_m_max <= 9.9:
                configuration['m_max_current'] = suggested_m_max
            id_do_again = False
    if suggested_m_max is None or suggested_m_max >= 9.9:
        print(f"!!!!! Procedure {METHODS[configuration['procedure_id']]} ({configuration['m_max_assessment']})")
        print(f"!!!!! could not asses correct m_max")
        print(f"!!!!! Final m_max is not set")
        print(f"!!!!! The current m_max value ({configuration['m_max_current']:4.2f}) remains")
        return
    # configuration['m_max'] = round(suggested_m_max, 2)
    configuration['m_max'] = round(configuration['m_max_current'], 2)
    configuration['sd_m_max'] = round(sd_m_max, 2)
    if 'lambda' in configuration and configuration.get('induced_seismicity', 'no') == 'yes':
        configuration['lambda_is'] = configuration['lambda'] * configuration['induced_seismicity_coefficient']
    compute_uncertainty(configuration)
    configuration['computation_time'] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def main():
    # warnings.filterwarnings("error", category=RuntimeWarning)
    configuration = load_configuration()
    # from ha3py.configuration import init_lambda_beta
    # init_lambda_beta(configuration, m_max=configuration['m_max_obs']+0.5)
    input_name = input_configuration_name()
    output_name = configuration.get('output_configuration')
    if input_name == output_name or output_name is None:
        if output_name is not None:
            print(f"Input and output configuration names '{output_name}' are the same")
        set_modify(configuration, 'output_configuration', 'Name of the output configuration file{} > ',
                                         " (other than '{0}' or enter to accept '{0}')", dtype=str)
    compute(configuration)
    print_separation_double_line()
    print('Final results')
    print_separation_single_line()
    names = get_events_occurrence(configuration).coefficient_names
    for name in names:
        print(f"{name:<10s} ={configuration[name]:7.3f} (+/- {configuration['sd_' + name]:6.3f} )")
    print_percentage_share(configuration)
    save_configuration(configuration)
    print_separation_double_line()
    print(f"Results written to the '{output_name}' file")
    print_separation_double_line()


if __name__ == "__main__":
    main()
