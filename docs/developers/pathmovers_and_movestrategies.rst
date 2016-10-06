.. _dev-pathmovers-movestrategies:

PathMovers and MoveStrategies
=============================

If you want to add new kinds of Monte Carlo moves, you need to know a little
more about the :class:`.PathMover` object. If you would like those Monte
Carlo moves to be compatible with the rest of the OPS path sampling
structure, you should learn more about :class:`MoveStrategies
<.MoveStrategy>`.


PathMovers
----------



EngineMovers
------------

Many of the new path sampling movers that people invent involve propagating
the dynamics in some way. To simplify this, we have an abstract object
called an :class:`.EngineMover`. Shooting moves are built based on the
:class:`.EngineMover`.

The big picture is that the :class:`.EngineMover` takes a :class:`.Sample`
from a given ``ensemble``, and propagates the trajectory can no longer
satisfy a given ``target_ensemble``. Details of subclasses depend mainly on
the ``ensemble``, the ``target_ensemble``, and the ``direction``. 

For example, shooting moves have the same ``ensemble`` and
``target_ensemble``. However, for extension moves (e.g., in the
:class:`.MinusMover`), the initial ``ensemble`` is different from the final
``target_ensemble``.

Before re-implementing anything, you should consider whether your needs are
met by the built-in subclasses of :class:`.EngineMover`. A few examples
might help to show how moves can be thought of as 

MoveStrategies
--------------

The :class:`.MoveStrategy` object acts as a factory for :class:`PathMovers
<.PathMover>`. This makes it much easier for the user to mix and match
various move types to create an overall custom scheme for the simulation.

