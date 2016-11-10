Custom Move Strategy and Path Mover
===================================

This example shows how to create a custom :class:`.MoveStrategy` and a
custom :class:`.PathMover`. Note that the custom path mover is very easy
here, but the custom move strategy is what makes it very easy to use with
the :class:`.MoveScheme` object, which facilitates analysis.

This particular example is on a simple toy model, and the new approach does
not seem to give much benefit. But this example also shows how to use tools
in OPS to compare, for example, replica travel time in two approaches.

-----

.. notebook:: examples/misc/custom_strategy_repex_shoot_repex.ipynb
   :skip_exceptions:
