.. _analysis-extractors:

Analysis Extractors
===================

Using extractors and extractor-filters
--------------------------------------

.. note::

    Extractor names are plural, even if they only extract a single item.
    This is for two reasons. First, the common use of extractors is expected
    to be to create an extractor-filter with the ``using`` method, in which
    case it reads more clearly in the plural (since it will iterate over
    many steps). Second, this allows a standard way to differentiate with
    filters. For example, the filter ``canonical_mover(mover)`` gives a
    filter that selects steps with the given mover/move type as the
    canonical mover. On the other hand, ``canonical_movers(step)`` returns
    the mover that is the canonical mover for that step.

Built-in extractors
-------------------

* ``trial_samples``
* ``active_samples``
* ``shooting_points``
* ``modified_shooting_points``
* ``canonical_movers``

* ``canonical_details``

Writing custom extractors
-------------------------

Extractors API
--------------

Extractor-Filters API
---------------------

In general, the way to create and use an extractor-filter is to create it
with the ``extractor.using(step_filter)`` method, and then to wrap the input
``steps`` with the resulting extractor-filter. Users will rarely directly
initialize an extractor-filter via its ``__init__`` method.
