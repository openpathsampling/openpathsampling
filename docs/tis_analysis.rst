.. _tis_analysis:

.. currentmodule:: openpathsampling.analysis.tis_analysis


TIS Analysis API
================

Abstract classes
----------------

.. autosummary::
   :toctree: api/generated/

    MultiEnsembleSamplingAnalyzer

Result data class
-----------------

.. autosummary::
   :toctree: api/generated/

    TransitionDictResults

Flux calculations
-----------------

.. autosummary::
   :toctree: api/generated/

    DictFlux
    MinusMoveFlux

Histogrammers
-------------

.. autosummary::
   :toctree: api/generated/

    EnsembleHistogrammer
    PathLengthHistogrammer

Transition probability
----------------------

.. autosummary::
   :toctree: api/generated/

    StandardTransitionProbability
    ConditionalTransitionProbability

Crossing probability function
-----------------------------

.. autosummary::
   :toctree: api/generated/
   
    FullHistogramMaxLambdas
    TotalCrossingProbability


Full TIS analysis
-----------------

.. autosummary::
   :toctree: api/generated/

    TISAnalysis
    StandardTISAnalysis
