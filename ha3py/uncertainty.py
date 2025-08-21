import numpy as np
from numpy.linalg import inv
from math import sqrt
from ha3py.ln_likelihood import ln_likelihood
from ha3py.get_events_occurrence import get_events_occurrence
from ha3py.constant_values import LN_10_
# from scipy.differentiate import hessian


# class FunctionForHessian:
#     def __init__(self, configuration):
#         self.configuration = configuration
#
#     def __call__(self, x_v):
#         out_np_arr = ln_likelihood(x_v, self.configuration)
#         return out_np_arr


def get_covariance(configuration, events_distribution=None, delta=0.01):
    if events_distribution is None:
        events_distribution = get_events_occurrence(configuration)
    x_v = events_distribution.coefficients[:-1]
    var_cov_size = len(x_v)
    dx_v = [x * delta for x in x_v]
    var_cov_inv_a = np.zeros((var_cov_size, var_cov_size))
    f_00 = ln_likelihood(x_v, configuration)
    for idx1 in range(var_cov_size):
        for idx2 in range(idx1, var_cov_size):
            x = x_v.copy()
            if idx1 == idx2:
                x[idx1] = x_v[idx1] - dx_v[idx1]
                f_bm1 = ln_likelihood(x, configuration)
                x[idx1] = x_v[idx1] + dx_v[idx1]
                f_bp1 = ln_likelihood(x, configuration)
                x[idx1] = x_v[idx1] - 2.0 * dx_v[idx1]
                f_bm2 = ln_likelihood(x, configuration)
                x[idx1] = x_v[idx1] + 2.0 * dx_v[idx1]
                f_bp2 = ln_likelihood(x, configuration)
                d1d1 = (-f_bp2 + 16 * f_bp1 - 30 * f_00 + 16 * f_bm1 - f_bm2) / (12 * dx_v[idx1] * dx_v[idx1])
                var_cov_inv_a[idx1][idx1] = - d1d1
            else:
                x[idx1] = x_v[idx1] - 1.0 * dx_v[idx1]
                x[idx2] = x_v[idx2] - 1.0 * dx_v[idx2]
                f_bmlm = ln_likelihood(x, configuration)
                x[idx1] = x_v[idx1] + 1.0 * dx_v[idx1]
                x[idx2] = x_v[idx2] - 1.0 * dx_v[idx2]
                f_bmlp = ln_likelihood(x, configuration)
                x[idx1] = x_v[idx1] - 1.0 * dx_v[idx1]
                x[idx2] = x_v[idx2] + 1.0 * dx_v[idx2]
                f_bplm = ln_likelihood(x, configuration)
                x[idx1] = x_v[idx1] + 1.0 * dx_v[idx1]
                x[idx2] = x_v[idx2] + 1.0 * dx_v[idx2]
                f_bplp = ln_likelihood(x, configuration)
                d1d2 = (f_bplp - f_bplm - f_bmlp + f_bmlm) / (4 * dx_v[idx1] * dx_v[idx2])
                var_cov_inv_a[idx1][idx2] = - d1d2
                var_cov_inv_a[idx2][idx1] = - d1d2
    var_cov_a = inv(var_cov_inv_a)
    for idx in range(var_cov_size):
        var_cov_a[idx][idx] = abs(var_cov_a[idx][idx])
    return var_cov_a



# def get_covariance(configuration, events_distribution=None, delta=0.01):
#     if events_distribution is None:
#         events_distribution = get_events_occurrence(configuration)
#     x_v = events_distribution.coefficients
#     if len(x_v) > 3:
#         return None
#     lamb = x_v[0]
#     beta = x_v[1]
#     d_beta = abs(beta) * delta  # d_beta   = 1% OF BETA   VALUE
#     d_lambda = abs(lamb) * delta  # d_lambda = 1% OF LAMBDA VALUE
#     f_00 = ln_likelihood(array([lamb, beta]), configuration)
#
#     # d2_Ln(LIKELIHOOD FUNCTION)/d2_BETA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#     f_bm2 = ln_likelihood(array([lamb, beta - 2 * d_beta]), configuration)
#     f_bm1 = ln_likelihood(array([lamb, beta - d_beta]), configuration)
#     f_bp1 = ln_likelihood(array([lamb, beta + d_beta]), configuration)
#     f_bp2 = ln_likelihood(array([lamb, beta + 2 * d_beta]), configuration)
#     d2d2_beta = (-f_bp2 + 16 * f_bp1 - 30 * f_00 + 16 * f_bm1 - f_bm2) / (12 * d_beta * d_beta)
#     # END OF d2_Ln(LIKELIHOOD FUNCTION)/d2_BETA DETERMINATION <<<<<<<<<<<<<<<<<
#
#     # d2_Ln(LIKELIHOOD FUNCTION)/d2_LAMBDA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#     f_lm2 = ln_likelihood(array([lamb - 2 * d_lambda, beta]), configuration)
#     f_lm1 = ln_likelihood(array([lamb - d_lambda, beta]), configuration)
#     f_lp1 = ln_likelihood(array([lamb + d_lambda, beta]), configuration)
#     f_lp2 = ln_likelihood(array([lamb + 2 * d_lambda, beta]), configuration)
#     d2d2_lambda = (-f_lp2 + 16 * f_lp1 - 30 * f_00 + 16 * f_lm1 - f_lm2) / (12 * d_lambda * d_lambda)
#     # END OF d2_Ln(LIKELIHOOD FUNCTION)/d2_LAMBDA DETERMINATION <<<<<<<<<<<<<<<
#
#     # d2_Ln(LIKELIHOOD FUNCTION)/d_BETA_LAMBDA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#     f_bmlm = ln_likelihood(array([lamb - d_lambda, beta - d_beta]), configuration)
#     f_bmlp = ln_likelihood(array([lamb + d_lambda, beta - d_beta]), configuration)
#     f_bplm = ln_likelihood(array([lamb - d_lambda, beta + d_beta]), configuration)
#     f_bplp = ln_likelihood(array([lamb + d_lambda, beta + d_beta]), configuration)
#     d2d_beta_lambda = (f_bplp - f_bplm - f_bmlp + f_bmlm) / (4 * d_beta * d_lambda)
#     # END OF d2_Ln(LIKELIHOOD FUNCTION)/d_BETA_LAMBDA DETERMINATION <<<<<<<<<<<
#
#     # VAR-COV MATRIX CALCULATION >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#     var_cov_inv_a = array([[-d2d2_lambda, -d2d_beta_lambda],
#                            [-d2d_beta_lambda, -d2d2_beta]])
#     var_cov_a = inv(var_cov_inv_a)
#     var_cov_a[0][0] = abs(var_cov_a[0][0])
#     var_cov_a[1][1] = abs(var_cov_a[1][1])
#     return var_cov_a


def compute_uncertainty(configuration, events_distribution=None):
    if events_distribution is None:
        events_distribution = get_events_occurrence(configuration)
    cov = get_covariance(configuration, events_distribution=events_distribution)
    # fun = FunctionForHessian(configuration)
    # x_v = events_distribution.coefficients[:-1]
    # cov1 = hessian(fun, x_v)
    if cov is None:
        return
    # cov1 = get_covariance1(configuration, events_distribution=events_distribution)
    configuration['cov_beta_lambda'] = cov.tolist()
    events_distribution_names = events_distribution.coefficient_names
    events_distribution_names = events_distribution_names[:-1]
    # configuration['coefficient_names'] = events_distribution_names
    for idx, name in enumerate(events_distribution_names):
        configuration[f'sd_{name}'] = sqrt(cov[idx][idx])
    if 'sd_beta' in configuration:
        configuration['sd_b'] = configuration['sd_beta'] / LN_10_
    if 'sd_lambda' in configuration\
            and configuration.get('induced_seismicity', 'no') == 'yes'\
            and 'induced_seismicity_coefficient' in configuration:
        configuration['sd_lambda_is'] = configuration['sd_lambda'] * configuration['induced_seismicity_coefficient']


if __name__ == "__main__":
    import sys
    from configuration import load_configuration, save_configuration
    params = load_configuration()
    ED = get_events_occurrence(params)
    compute_uncertainty(params, events_distribution=ED)
    ED_names = ED.coefficient_names
    for parameter_name in ED_names:
        sd_name = 'sd_' + parameter_name
        if sd_name not in params:
            print(f'Missing {sd_name})')
        elif parameter_name not in params:
            print(f'Missing {parameter_name})')
        else:
            print(f'{parameter_name}: {params[parameter_name]} (+/- {params[sd_name]})')
    if len(sys.argv) >= 3:
        save_configuration(params, sys.argv[2])
