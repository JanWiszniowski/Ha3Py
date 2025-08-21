.. _configuration:

#########################
Configuration description
#########################

The configuration is kept in the Python dictionary, where keys are strings and values are depend on parameters:
strings, float values, integer values, boolean values, sub-dictionaries, or lists.
The configuration file (example name: ``Ha3Py.json``) is a file in JavaScript Object Notation (JSON).
Here is the example file::

    {
        "paleo_catalog": {
            "begin": -10000.0,
            "end": 1800.0,
            "time_span": 11800.0,
            "sd": 0.5,
            "m_min": 5.94,
            "name": "p.txt",
            "earthquakes": [
                {
                    "magnitude": 6.02,
                    "date": -9650.0,
                    "time_span": 350.0,
                    "sd": 0.5
                },
                {
                    "magnitude": 6.01,
                    "date": -8887.0,
                    "time_span": 763.0,
                    "sd": 0.25
                },
    ...
    ...
                {
                    "magnitude": 6.08,
                    "date": 1203.0,
                    "time_span": 247.0,
                    "sd": 0.25
                }
            ]
        },
        "historic_catalog": {
            "begin": 1550.0,
            "end": 1800.0,
            "time_span": 250.0,
            "m_min": 3.5,
            "name": "h.txt",
            "earthquakes": [
                {
                    "magnitude": 4.16,
                    "date": 1551.704109589041,
                    "time_span": 1.7041095890410816,
                    "sd": 0.25
                },
                {
                    "magnitude": 4.97,
                    "date": 1553.7452054794521,
                    "time_span": 2.041095890411043,
                    "sd": 0.25
                },
    ...
    ...
                {
                    "magnitude": 4.59,
                    "date": 1799.5205479452054,
                    "time_span": 3.263013698630175,
                    "sd": 0.25
                }
            ]
        },
        "complete_catalogs": [
            {
                "begin": 1800.0,
                "end": 1900.0,
                "time_span": 100.0,
                "m_min": 3.5,
                "sd": 0.1,
                "name": "c1.txt",
                "earthquakes": [
                    {
                        "magnitude": 3.61,
                        "date": 1850.0,
                        "sd": 0.1
                    },
                    {
                        "magnitude": 4.86,
                        "date": 1850.0,
                        "sd": 0.1
                    },
    ...
    ...
                    {
                        "magnitude": 3.03,
                        "date": 1950.0,
                        "sd": 0.1
                    }
                ]
            }
        ],
        "begin": -10000.0,
        "end": 2000.0,
        "output_text_file": "Ha3PyMTest",
        "area_name": "Area",
        "created_on": "2023-06-04 09:42:14",
        "m_max_obs": 6.08,
        "sd_m_max_obs": 0.25,
        "induced_seismicity": "no",
        "induced_seismicity_coefficient": 1.0,
        "time_intervals": [
            1.0,
            50.0,
            470.0,
            1000.0
        ],
        "time_span": 12000.0,
        "m_min": 3.0,
        "procedure id": 8,
        "prior_m_max": 6.08,
        "sd_prior_m_max": 0.5,
        "m_max_prior": 7.0,
        "sd_m_max_prior": 0.25,
        "prior_b": "",
        "prior_beta": null,
        "sd_prior_beta": null,
        "likelihood_optimization_method": null,
        "m_max_current": 8.99,
        "q_beta": 16.0,
        "q_lambda": 16.0,
        "beta": 2.619477443131598,
        "lambda": 7.878695391805496,
        "magnitude_distribution": "GutenbergRichterBayes",
        "occurrence_probability": "Poisson",
        "delta": "Kijko-Sellevoll",
        "m_max_suggested": 7,
        "sd_m_max": 0.5,
        "cov_beta_lambda": [
            [
                0.0041074392844434576,
                -0.03865960610149208
            ],
            [
                -0.03865960610149208,
                1.050816279158931
            ]
        ],
        "sd_beta": 0.06408930709910553,
        "sd_lambda": 1.0250933026602658,
        "b": 1.0516286240033526,
        "sd_b": 0.02783363242214443,
        "lambda_is": 9.532098273353775,
        "sd_lambda_is": 1.0250933026602658,
        "hazard": [
            {
                "mag": 3.0,
                "lambda_mag": 7.703958463289057,
                "return_period": 0.12980339974121163,
                "probabilities": [
                    1.0,
                    1.0,
                    1.0,
                    1.0
                ]
            },
            {
                "mag": 3.1,
                "lambda_mag": 5.8736839448960225,
                "return_period": 0.1702509037567397,
                "probabilities": [
                    0.9932831012788538,
                    1.0,
                    1.0,
                    1.0
                ]
            },
    ...
    ...
            {
                "mag": 8.900000000000006,
                "lambda_mag": 1.2599577011667585e-05,
                "return_period": 79367.74378012613,
                "probabilities": [
                    1.259949267673477e-05,
                    0.0006297680614141266,
                    0.005903212762204291,
                    0.012515638410797392
                ]
            }
        ],
        "m_max": 8.99
    }

Parameters description
######################

Below is description of parameters (not all described bellow parameters are in the example)

Catalogues
==========

Three types of catalogues can be defined. The paleo-catalogue and historical catalogue
can be defined only once, whereas many complete catalogues having different period and
completeness magnitude can be defined. Catalogues must be declustered
and can not contain event with magnitude below the defined minimum magnitude :math:`m_{min}`.
Magnitudes in all catalogs must be unified.

:paleo_catalog: (dict) The prehistoric (paleo-) catalogue
:historic_catalog: (dict) The historic catalogue
:complete_catalogs: (list(dict)) The list of complete (instrumental) catalogues

Catalogue parameters
--------------------

:begin: (float) The begin of the catalogue time in years AD.
:end: (float) The end of the catalogue time in years AD.
:time_span: (float) The catalogue time span in years.
:sd: (float) The standard magnitude deviation of earthquakes in the catalogue.
    This value is valid unless individual sd is defined for specific earthquake.
:m_min: (float, str) The minimum (completeness) magnitude :math:`m_{min}`.
    If it is float value, it means the :math:`m_{min}`.
    WARNING. The catalogues can not contain earthquakes with magnitude below :math:`m_{min}`.
    If *m_min* is the text '*variable*',
    it mean that the completeness level is variable and is defined with events.
    The '*variable*' option is allowed only to complete catalogues.
:name: (str) The name of the catalogue. Any name is allowed, though it depends on the imported catalogue format.
:earthquakes: (list(dict)) The list o earthquakes parameters described bellow.

Earthquake parameters
---------------------

Specific earthquake parameters depends on the catalogue type.

:magnitude: (float) The earthquake magnitude. Magnitudes in all catalogs must be unified.
:date: (float) The date of the earthquake in years AD.
:time_span: (float) The time_span to the previous event.
    It is required in historical and prehistorical catalogues.
:sd: (float) The standard deviation of the earthquake magnitude.
:m_min: (float) The minimum catalogue magnitude attributable to the period
    of occurrence of this earthquake.
    It is required in complete catalogue when *m_min* of the catalogue is '*variable*'.

Earthquake occurrence and magnitude distribution parameters and coefficients
============================================================================

Earthquake occurrence and magnitude distribution parameters
-----------------------------------------------------------

:magnitude_distribution: (str) The name of magnitude distribution model e.g 'GCompound Gutenberg-Richter',
:occurrence_probability: (str) The name of earthquake occurrence model e.g 'Poisson-gamma compound',

Earthquake occurrence and magnitude distribution coefficients
-------------------------------------------------------------

Earthquake occurrence and magnitude distribution coefficients depend on the
earthquake occurrence and magnitude distribution classes used for seismic hazard assessment
and defined by the parameters.
The example below refers to a typical solution::

    "beta": 2.619477443131598,
    "lambda": 7.878695391805496,
    "q_beta": 16.0,
    "q_lambda": 16.0,
    "prior_beta": null,
    "sd_prior_beta": null,
    "sd_beta": 0.06408930709910553,
    "sd_lambda": 1.0250933026602658,

``prior_*`` (where ``*`` means the coefficient name) coefficient are required,
if we use bayesian approach to the earthquake occurrence assessment.

Maximum magnitude parameters
============================

:procedure_id: (int) The number of the selected in the ``ha_config`` or ``ha3`` programs maximum assessment method.
:m_max_method: (str) The name of the defined maximum assessment method.
    It is set in the ``ha_config`` or ``ha3`` programs and used on for report.
    If you configure manually, you should define or redefine it to have the correct report.
:m_max_obs: (float) The maximum ever observed in the studied region magnitude.
:sd_m_max_obs: (float) The standard deviation of the maximum observed magnitude.
:prior_m_max: (float) The prior maximum magnitude used in bayesian maximum magnitude assessments.
:sd_prior_m_max: (float) The standard deviation of the prior maximum magnitude
    used in bayesian maximum magnitude assessments.
:m_max_suggested: (float) The temporary suggested maximum magnitude
    which is the result od maximum magnitude computation.
    It is not final maximum magnitude.
:m_max_current: (float) The temporary suggested maximum magnitude used for
    magnitude distribution definition used for computation before final maximum magnitude is set.
    When the final maximum magnitude is set, the m_max_current should equals the final maximum magnitude.
:m_max: (float) The final maximum magnitude.
:sd_m_max: (float) The standard deviation of the final maximum magnitude.
:m_max_assessment: (str) The non bayesian maximum magnitude assessment method.
:bayesian_m_max_assessment: (str) The bayesian maximum magnitude assessment method.
    When this parameter is not defined,
    the maximum magnitude is estimated by the method defined in the 'm_max_assessment' parameter.
    When this parameter is not defined,
    the maximum magnitude is estimated by this method
    whereas the 'm_max_assessment' parameter is used in some bayesian methods,
    which requires the pre-estimation of the maximum magnitude from data.
:bayesian_m_max_estimator: (str) The bayesian maximum magnitude likelihood point estimator.
    Three estimators are available: 'max' - maximum likelihood, 'expected' - expected likelihood,
    and 'median' - median likelihood treated as a probability distribution.
:m_min_ref: (float) The minimum magnitude definition defined for maximum magnitude estimation.
    It is optional. When is not defined, the minimum magnitude of all catalogs is used ``m_min``.
    However, ``m_min`` can be small, which is numerically not recommended.

Other parameters
================

:begin: (float) The begin of all catalogues time in years AD.
:end: (float) The end of all catalogues time in years AD.
:time_span: (float) The all catalogues time span in years.
:output_text_file: (str) The name of the output text with the estimation results report.
:area_name: (str) The name of the investigated area
:created_on: (str) The date and time of data creation e.g."2023-06-04 09:42:14",
:induced_seismicity: (str) Information if the investigated seismicity
    is the case of induced seismicity ("no" or "yes")
:induced_seismicity_coefficient: (float),
:time_intervals: (list)
:m_min: (float),
:prior_b: "",
:likelihood_optimization_method: (str) The optimisation method used for maximum likelihood computation.
:COV: (np.array) The covariance matrix (Inversion of ln ln_likelihood hestian).
:coefficients: (list(str)) List if earthquake occurrence probability estimated coefficients names,
    e.g. ['lambda', 'beta', 'm_max']
    (see :ref:`Earthquake occurrence and magnitude distribution coefficients`).
:b: (float) The b-value of Gutenberg-Richter magnitude distribution,
:sd_b: (float) The standard deviation of the b-value,
:lambda_is: (float) The annual :math:`\lambda` in the case of induced seismicity,
:sd_lambda_is: (float) The standard deviation of the induced seismicity :math:`\lambda`,
:hazard: (list(dict)) The list of magnitude exceedance probabilities, return periods, etc.
    for various magnitudes (see :ref:`Seismic hazard values`).

Seismic hazard values
=====================

:mag: (float) The magnitude value, for which the following parameters were calculated.
:lambda_mag: (float)  ???? 1.2599577011667585e-05,
:return_period: (float) the return period in years
:probabilities: (list(float)) List of non exceeding of above magnitude *mag* probability
    in periods defined in the *time_intervals* list

Simulation configuration
========================

:simulation: (dict)

Simulation parameter are analogues to catalogs parameters.
The dictionary `simulation` consists of:

:pre-historic data: (dict) which has description, how the generate the pre-historic catalogue,
:historic data: (dict) which has description, how the generate the historic catalogue,
:complete data: (list) which is the list of a few complete catalogues,
    each has the description, how the generate the complete catalogue.

Description, how the generate the catalogue consists of:

:generator: (str) The generator name. There are allowed three names:
    'full_simulation_incremental', 'extreme_simulation', and 'no_date_simulation'.
    For paleo- and historic catalogues 'full_simulation_incremental' and 'extreme_simulation'
    can be applied, wheras for complete catalogues 'extreme_simulation'
    and 'no_date_simulation' can be applied.
:begin: (float) The begin of the catalogues in years AD.
:end: (float) The end of the catalogues in years AD.
:time_span: (float) The time interval of the catalogue
:time_uncertainty: The event time uncertainty in the catalogue
:magnitude_uncertainty: (float) The magnitude uncertainty in the catalogue
:name: (str) The name of the catalogue
:m_min: (float) The minimum magnitude in the catalogue (required)
