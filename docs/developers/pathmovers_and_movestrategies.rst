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
might help to show how moves can be thought of as including an
:class:`.EngineMover`, rather than needing to subclass one:

TODO

MoveStrategies
--------------

The :class:`.MoveStrategy` object acts as a factory for :class:`PathMovers
<.PathMover>`. This makes it much easier for the user to mix and match
various move types to create an overall custom scheme for the simulation.

When creating a new move strategy, you'll mainly need to override the
:meth:`make_movers() <MoveStrategy.make_movers>` method. You'll also need to
decide the ``level`` of the strategy. The different levels are essentially
priorities. When a :class:`.MoveScheme` builds its movers, the strategies
are applied in an order sorted first by ``level``, and second by order the
strategy was added to the scheme. Identifying the correct ``level`` for you
strategy is by far the most complicated part of adding a move strategy, so
the rest of this section will explain it.

Let's start with the way that movers are organized with a
:class:`.MoveScheme`. Every scheme has organizes its movers into "groups."
These groups should define the "canonical" move types (that is, the way you
normally think about a move in path space: shooting, replica exchange, etc.)
The process of building a move scheme goes in the following order:

1. For each group, decide what ensembles will be involved as input and
   output (``SIGNATURE`` level).
2. Create the specific path movers for those ensembles (``MOVER`` level).
3. Create the groups that will organize the move scheme. (``GROUP`` level).
4. If necessary, change the kind of mover or reorganize the group
   (``SUPERGROUP`` level).
5. Finally, create the global organization of the move scheme (organize the
   groups themselves, ``GLOBAL`` level).

The differences between ``GROUP`` and ``SUPERGROUP`` are pretty flexible,
and many strategies could work for either one. Here are some examples of
when to use each level:

* ``levels.SIGNATURE``: Use this if your mover changes *which ensembles* are
  involved. For example, different replica exchange strategies include
  "nearest neighbor," "all possible," and "specific selected." These don't
  change the nature of the move, but do change which ensembles are involved.
* ``levels.MOVER``: Use this if the nature of the move changes, but not its
  fundamental purpose. For example, if you are implemented a different kind
  of shooting strategy (a different approach for shooting point selection,
  or two-way shooting instead of one-way), this would be the correct
  approach. It changes the movers without changing the input and output
  ensembles of each move.
* ``levels.GROUP``: Use this if you're creating a new group of movers from
  movers in other groups. For example, rather than randomly selected which
  moves to do from a group, you might create a new group where you combine
  them 
* ``levels.SUPERGROUP``: Use this if you're 
* ``levels.GLOBAL``: Use this to organize the global structure of 

In most cases, you'll probably be adding a new type of mover. In that case,
you should use ``levels.MOVER``. The best approach is, as always, to find an
example in the code that does something similar to what you want to do.

Technically, each level is associated with an integer value, and you can add
other levels between (much like Python's ``logging`` facilities). The
``levels`` object just gives convenient access to specific values (10, 30,
50, 70, 90). However, we don't recommend straying from those default levels
unless you're very certain that you must.


