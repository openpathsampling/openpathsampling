.. _analysis-filters:

Analysis Filters
================

For a detailed overview of the using extractors and filters to simplify
analysis functions, see ???. This page will document built-in filters, as
well as the underlying API for creating new filters.


Built-in filters
----------------

Step filters
~~~~~~~~~~~~

* ``all_steps``
* ``accepted_steps``
* ``rejected_steps``
* ``shooting_steps``

* ``canonical_mover``
* ``trial_replica``
* ``trial_ensemble``


Sample filters
~~~~~~~~~~~~~~

Frequently, you'll extract a list of samples from the step, e.g., using
either the ``active_samples`` extractor or the ``trial_samples`` extractor.
In this case, you may want to select only specific samples (perhaps based on
the ensemble) for your analysis. Sample filters make this possible, and can
be used as the ``secondary_filter`` when creating extractor-filters.

* ``all_samples``

* ``ensemble``
* ``replica``
* ``sampling_ensemble``
* ``minus_ensemble``
* ``ms_outer_ensemble`


Writing custom filters
----------------------

Filters API
-----------



