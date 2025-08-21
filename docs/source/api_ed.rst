.. _api_ed:

Seismic event occurrence probability
====================================

The seismic event occurrence probability define the probability of
occurrence of N events in the defined time duration.
All algorithms work on an abstract event occurrence probability class `BaseNEventsOccurrenceProb`.
It allows the assessment of various event occurrence probability,
which depend on various parameters.

The special abstract subclass of the `BaseNEventsOccurrenceProb` is the `LambdaNEventsOccurrence`,
which should be the base class of all event occurrence probability classes
that have only one :math:`\lambda` variable parameter.

Base event occurrence probability classes
=========================================

.. automodule:: BaseOccurrence
   :members:

.. automodule:: LambdaOccurrence
   :members:

Users can use one from two event occurrence probability classes:
* classic Poisson probability,
* gamma compound Poisson probability.

Predefined events event occurrence probability classes
======================================================

.. automodule:: PoissonOccurrence
   :members:

.. automodule:: GammaPoissonOccurrence
   :members:

