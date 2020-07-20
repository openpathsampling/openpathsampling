.. _network:


Networks and Transitions
========================

Networks are an easy way to set up ensembles for path sampling calculations,
and to retain a context for those ensembles when analyzing the results of a
simulation. Networks are made up of multiple transitions, where each
transition represents a directed :math:`A\to B` reaction.

Many path sampling methods (transition interface sampling as one example)
require sampling multiple path ensembles, and then results from those path
ensembles are combined in analysis. As such, the path ensemble itself is not
enough: you need both the path ensemble and its context. Networks provide
that context.

Networks are so named because they allow the study of complicated transition
networks, instead of just single transitions. It is at this level that we
encounted ideas like multiple state path sampling methods and multiple
interface set transition interface sampling.

A longer discussion of networks, and the associated concept of transitions,
can be found in the file `Transitions and Networks
<guides/transitions_and_networks.html>`_.

.. currentmodule:: openpathsampling.high_level.network

Abstract network class
----------------------
.. autosummary::
   :toctree: api/generated

   TransitionNetwork

TPS networks
------------
.. autosummary::
   :toctree: api/generated

   TPSNetwork
   FixedLengthTPSNetwork

TIS networks
------------
.. autosummary::
   :toctree: api/generated

   MSTISNetwork
   MISTISNetwork

-----

Transitions
-----------

.. currentmodule:: openpathsampling.high_level.transition

In general, users don't need to directly create a :class:`.Transition` --
the network acts as a factory that creates them.

.. autosummary::
   :toctree: api/generated

   Transition
   TPSTransition
   FixedLengthTPSTransition
   TISTransition

