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

The non-parametric magnitude distribution is sensitive to the completeness of the magnitudes.
Therefore, applying this to complete and extreme catalogues should be done with caution,
as additional conditions must be fulfill, e.g., there must be no gaps in the magnitude ranges.

Predefined magnitude distribution classes
-----------------------------------------

.. automodule:: GutenbergRichter
   :members:

.. automodule:: CompoundGutenbergRichter
   :members:

.. automodule:: NonParametricGaussian
   :members:
