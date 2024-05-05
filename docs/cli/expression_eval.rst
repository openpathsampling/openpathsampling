.. _expression_eval:

Expression Evaluation
=====================

The ``compile`` command has support for evaluation of custom expressions.

In general, these tools are a very limited subset of Python, but they do
make some this a little easier. For example, if you want the edge of a state
definition to be at 100 degrees, but the underlying CV actually calculates
in radians, you can use the expression value::

    100 * np.pi / 180

Internally, OpenPathSampling will evaluate this and return 1.745329(...),
which is 100 degrees in radians.

Expression evaluation is limited to a single `Python expression
<https://docs.python.org/3/reference/expressions.html>`_. Typically, this
translates to "a single line of Python." No imports are allowed, but the
following namespaces are available:

* ``np``: NumPy, as if ``import numpy as np``
* ``math``: standard math library, as if ``import math``

.. _EvalFloat:

EvalFloat
---------

A parameter of type ``EvalFloat`` uses expression evaluation and must
evaluate to a floating point number. It will be cast to a Python ``float``,
and if that cast fails, an error will be raised.

.. _EvalInt:

EvalInt
-------

A parameter of type ``EvalInt`` uses expression evaluation and will be cast
to a Python ``int`` (i.e., the result will be passed through
``int(result)``). Note that this means that it will use the standard Python
truncation rules for ``int`` (``int(1.5) == 1``, ``int(-1.5) == -1``, etc.).
If that cast fails, an error will be raised.

.. _EvalIntStrictPos:

EvalIntStrictPos
----------------

A parameter of type ``EvalIntStrictPos`` follows the same rules as a
parameter of type :ref:`EvalInt <EvalInt>` with the additional restriction
that the resulting value must be strictly positive (non-negative, non-zero).
If it is not, an error will be raised.
