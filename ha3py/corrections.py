from math import sqrt, exp


def add_catalog_var(catalogue, sum_var, nr_eq):
    catalogue_sd = catalogue.get('sd')
    for eq_phs in catalogue['earthquakes']:
        sum_var += eq_phs.get('sd', catalogue_sd) ** 2
        nr_eq += 1
    return sum_var, nr_eq


def lambda_correction(configuration):
    r"""
    Required parameters:
        paleo_catalog
        historic_catalog
        complete_catalogs
    """
    lamb = configuration.get('lambda')
    if lamb is None:
        return
    beta = configuration.get('beta')
    if beta is None:
        return
    p_phs = configuration.get('paleo_catalog')
    p_his = configuration.get('historic_catalog')
    p_comp = configuration.get('complete_catalogs')
    m_min = configuration['m_min']
    m_max = configuration['m_max_current']
    x = m_min
    nr_eq = 0
    sum_var = 0.0
    # ======================================================================
    if p_phs:  # CONTRIBUTION TO VAR(MAGNITUDE) FROM PRE-HISTORIC DATA
        sum_var, nr_eq = add_catalog_var(p_phs, sum_var, nr_eq)
    if p_his:  # CALCULATION OF CONTRIBUTION TO VAR(MAG) FROM HISTORIC DATA
        sum_var, nr_eq = add_catalog_var(p_his, sum_var, nr_eq)
    if p_comp:  # CALCULATION OF CONTRIBUTION TO VAR(MAG) FROM COMPLETE DATA
        for complete_catalog in p_comp:
            sum_var, nr_eq = add_catalog_var(complete_catalog, sum_var, nr_eq)
    # ======================================================================
    # CORRECTION FACTOR CALCULATION (eq. 21), A.Kijko & M.A. Sellevoll, 1992,
    # "Estimation of earthquake hazard parameters from incomplete data files,
    # Part II. Incorporation of magnitude heterogeneity", Bull. Seism. Soc.
    # Am., Vol.82, pp 120-134.
    var_mag = sum_var / nr_eq  # THE MEAN VARIANCE OF EARTHQUKE
    # MAGNITUDE DETERMINATION
    sd_mag = sqrt(var_mag)  # THE MEAN SD OF EARTHQUKE
    # MAGNITUDE DETERMINATION
    xx = x
    if xx < m_min:
        xx = m_min
    if xx > m_max:
        xx = m_max
    sqrt_2_ = sqrt(2)
    gamma = beta * sd_mag / sqrt_2_
    gamma2 = gamma ** 2
    # cf = 0.5 * exp(-gamma2) * (1+erf((mag_max-xx)/(1.4142*sd_mag)+gamma)); ?????
    # return 0.5 * exp(gamma2) * (1 + erf((mag_max-xx)/(sqrt_2_*sd_mag)+gamma))
    configuration['lambda'] = lamb * exp(-gamma2)
# END OF CORRECTION FACTOR CALCULATION -------------------------------
# ======================================================================
