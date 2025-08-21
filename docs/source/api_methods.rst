.. _api_methods:

Seismic event occurrence and magnitude assessment methods
#########################################################

Catalogues import and program configuration methods
===================================================

.. automodule:: configuration
   :members:

Creating Ha3Py objects or selecting methods based on their names
================================================================

.. automodule:: get_events_occurrence
   :members:

.. automodule:: get_magnitude_distribution
   :members:

`---------------------------------------------------------------------------`

The creating :math:`\Delta` computation objects method is described with the
:ref:`Della classes <api_delta>` description.

Magnitude occurrence probability assessment methods
===================================================

.. automodule:: compute
   :members:

.. automodule:: ln_likelihood
   :members:

The :math:`m_{max}` estimation procedures
=========================================

These procedures estimate the maximum possible magnitude in the area.
They all use :math:`m_{max}^{obs}`,
but some procedures require the earthquake occurrence probability,
delta estimation methods, or catalogues.

.. automodule:: fsolve_delta
   :members:

.. automodule:: iteration
   :members:

.. automodule:: momentum
   :members:

.. automodule:: gibowicz_kijko
   :members:

.. automodule:: robson_whitlock
   :members:

.. automodule:: bayesian_normal
   :members:

.. automodule:: bayesian_by_shift
   :members:

.. automodule:: bayesian_fiducial
   :members:

.. automodule:: primitive
   :members:

The :math:`m_{max}` support procedures
======================================

.. automodule:: compute_m_max
   :members:

.. automodule:: m_max_utils
   :members:

Uncertainty assessment methods
==============================

.. automodule:: uncertainty
   :members:

Correction procedures
=====================

.. automodule:: corrections
   :members:

Result visualization
====================

.. automodule:: return_period
   :members:

.. automodule:: print_results
   :members:

.. automodule:: plot_results
   :members:


