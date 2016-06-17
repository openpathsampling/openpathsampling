.. _common-setups:

=============
Common Setups
=============

OpenPathSampling is very flexible and simplifies many different approaches
to path sampling. Since it is easy to get lost under all those options, here
we summarize the most common setups of networks and move schemes for
transition path sampling and transition interface sampling.


Throughout the following, we assume that ``A``, ``B``, and ``C`` are
volumes that represent states. We assume that ``interfacesA`` is a list of
interfaces leaving state ``A``, and that ``interfacesAB`` is a list of
interfaces leaving state ``A`` in the expected direction of ``B``.
Similarly, we expect ``orderparameterA`` is a collective variable describing
the process of leaving ``A``, and ``orderparameterAB`` is a CV for leaving
``A`` toward ``B``.

------------------------
Transition Path Sampling
------------------------

Transition networks for TPS
===========================

For all transition path sampling setups, you can use either fixed path
length ensembles or flexible path length ensembles. The fixed path length
ensembles are created with the :class:`.FixedLengthTPSNetwork`, while the
flexible path length ensembles are created with the :class:`.TPSNetwork`.
The interfaces to these classes are the same, except that any method to
create a :class:`.FixedLengthTPSNetwork` requires an extra ``length``
parameter.  The same network classes are used for both multiple state TPS
and 2-state TPS.

Unidirectional two-state TPS
----------------------------

The standard initialization of the TPS network objects takes the parameters
``initial_states`` and ``final_states``. If you want a :class:`.TPSNetwork`
from state ``A`` to state ``B``, that can be created with:

.. code-block:: python

   network = paths.TPSNetwork(initial_states=A, final_states=B)

This will only sample the :math:`A\to B` transitions. The TPS examples for
alanine dipeptide give practical illustrations of this for both fixed-length
and flexible-length ensembles.

Bidirectional two-state TPS
---------------------------

If you want to sample both the :math:`A\to B` transitions and the
:math:`B\to A` transitions in one TPS simulation, then you can achieve the
desired sampling ensemble in several ways.

First, you could use the approach as above, but using
``initial_states=[A,B]`` and ``final_states=[A,B]``. But there is a shortcut
for this (particularly useful when there are many states) with

.. code-block:: python

   network = paths.TPSNetwork.from_states_all_to_all([A,B])

Both of these approaches ignore the self-transitions (:math:`A\to A` and
:math:`B\to B`). You can also explicitly create the transitions with the
:meth:`.TPSNetwork.from_state_pairs` method:

.. code-block:: python

   network = paths.TPSNetwork.from_state_pairs([(A,B), (B,A)])

Multiple-state TPS
------------------

The approaches for bidirectional 2-state TPS can be directly generalized to
multiple states. By default, the multiple state ensemble ignores
self-transitions (:math:`A\to A`) if created with either the standard
initialization or with :meth:`from_states_all_to_all
<.TPSNetwork.from_states_all_to_all>`.  If you would
like to include the self-transitions (not likely), then add the parameter
``allow_self_transitions=True`` to either the standard initialization or to
``from_states_all_to_all``. 

The method :meth:`.from_state_pairs` *will* allow self-transitions if you
explicitly include them. In practice, :meth:`.from_state_pairs` is most
useful when you only want to include certain transitions. For example, if
you have states ``A``, ``B``, and ``C``, but only want to include the
transitions :math:`A\to B`, :math:`A\to C`, and :math:`B\to C`, this is the
method to use. You'd create the network with

.. code-block:: python

   network = paths.TPSNetwork.from_state_pairs([(A,B), (A,C), (B,C)])

Move schemes for TPS
====================

Often, when using TPS (and especially when using flexible-length TPS), the
entire move scheme consists of a single shooting mover. Currently, OPS only
supports one-way shooting. A move scheme consisting of a single one-way
shooting move can be created with the :class:`.OneWayShootingMoveScheme`.
The common way to set this up is:

.. code-block:: python

   scheme = paths.OneWayShootingScheme(network, selector, engine)

where ``network`` would be the (TPS) network, ``selector`` is a shooting point
selector (usually an instance of :class:`.UniformSelector`), and ``engine``
is the desired dynamics engine.

*****

-----------------------------
Transition Interface Sampling
-----------------------------

Transition networks for TIS
===========================

As with TPS, the 2-state system in TIS is just a special case of the
multiple-state approach, so we use the same network classes to create
2-state systems as multiple-state systems.

.. _unidirectional-TIS:

Unidirectional two-state TIS
----------------------------

For unidirectional 2-state TIS (only studying the transition :math:`A\to B`,
not :math:`B\to A`), we use the :class:`.MISTISNetwork` with only one
transition listed:

.. code-block:: python

   network = paths.MISTISNetwork([(A, interfacesAB, orderparameterAB, B)])

This will sample the transition from ``A`` to ``B`` using the list of
``interfaces``, and the resulting analysis will be based on the collective
variable ``orderparameter``.

.. _bidirectional-TIS:

Bidirectional two-state TIS
---------------------------

For bidirectional 2-state TIS (simultaneously studying both the :math:`A\to
B` transition and the :math:`B\to A` transition), you could use a
:class:`.MISTISNetwork` as in the unidirectional case, but giving both
transitions instead. However, using the :class:`.MSTISNetwork` is a little
simpler, and gives completely equivalent results:

.. code-block:: python

   network = paths.MSTISNetwork([(A, interfacesA, orderparameterA),
                                 (B, interfacesB, orderparameterB)])


Multiple-state TIS
------------------

The network for the standard multiple state TIS, where there is one set of
interfaces for each state, is given by straightforward extension of the
bidirectional 2-state case to more states. Illustrations of this are in the
toy model MSTIS example, and in the alanine dipeptide MSTIS example.

If you wish to only focus on certain final states from a particular initial
states, or if you want to use more than one interface set per initial state,
then you need to use the multiple interface set variant of multiple state
TIS. This is given by straightforward extension of the unidirectional case
to more transitions. An illustration of this is in the toy model MISTIS
example.

More details on the distinctions between the MSTIS network and the MISTIS
networks are in the :ref:`"Which network should I use" <which-network>`
section.

Move schemes for TIS
====================

Here we'll discuss some standard and simple move schemes for TIS, which tend
to be significantly more complicated than for TPS. If you want a much more
complicated move scheme, it is usually good to start with one of these basic
move schemes, and then to use the :class:`.MoveStrategy` objects to modify
the scheme.

Standard TIS scheme
-------------------

The default TIS scheme includes one-way shooting (uniform shooting point
selection) for each TIS and multiple state ensemble, path reversal movers on
those same ensembles, a minus mover for each state, and nearest-neighbor
replica exchange. The probabilities of choosing each move type are designed
such that, for each ensemble, path reversal and replica exchange are tried
half as frequently as shooting. The minus move is tried 1/5 as frequently as
shooting.

This move scheme is generated with

.. code-block:: python

   scheme = paths.DefaultScheme(network, engine)

Single replica TIS
------------------

Any move scheme can be converted to a single replica move scheme with ???
