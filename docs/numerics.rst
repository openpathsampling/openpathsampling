.. _numerics:

.. currentmodule:: openpathsampling.numerics

Numerics
========

While most numerical functions needed by OPS are provided by libraries such
as numpy, there are a few specialized tools we have implemented. These are
in the numerics subpackage.


Resampling
----------

Tools for getting errors on pandas DataFrames.

.. autosummary::
   :toctree: api/generated/

    ResamplingStatistics
    BlockResampling


Lookup Functions
----------------

Interpolation tools that turn tables into functions.

.. autosummary::
   :toctree: api/generated/

    LookupFunction
    LookupFunctionGroup
    VoxelLookupFunction

Histograms
----------

.. autosummary::
   :toctree: api/generated/

    Histogram
    SparseHistogram
    HistogramPlotter2D
    histograms_to_pandas_dataframe
    Histogrammer


Histogram Combiners
-------------------

.. autosummary::
   :toctree: api/generated/

    WHAM

