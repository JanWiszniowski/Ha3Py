.. _api_md:

Magnitude distribution
======================

The magnitude distribution define the probability of exceeding the magnitude.
The magnitude distribution realises an object-oriented approach to this issue.
All algorithms work on an abstract magnitude probability class `MagnitudeDistribution`.
It allows the assessment of various magnitude distribution,
which depend on various parameters.

Base magnitude distribution class
---------------------------------

.. automodule:: MagnitudeDistribution
   :members:

Users can use one from three predefined magnitude distribution classes:
* classic Gutenberg-Richter magnitude distribution,
* compound Gutenberg-Richter magnitude distribution,
* non-parametric magnitude distribution.

Predefined magnitude distribution classes
-----------------------------------------

.. automodule:: GutenbergRichter
   :members:

.. automodule:: CompoundGutenbergRichter
   :members:

.. automodule:: NonParametricGaussian
   :members:
