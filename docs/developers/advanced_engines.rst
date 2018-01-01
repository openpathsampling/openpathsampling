.. _advanced_engines:

Advanced Engines
================

Custom snapshot features
------------------------

We frequently point to the fact that OPS can be used for any kind of
dynamics, even things that have nothing to do with molecular dynamics. If
that's your interest (or if you're interested in methods like *ab initio* or
quantum dynamics which might make use of a definition of state beyond
standard MD), then this is the section for you.

Standard snapshot features
~~~~~~~~~~~~~~~~~~~~~~~~~~

Adding features with new stored data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Adding features that are proxies to stored data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Adding features that are properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Support for additional OPS modules
----------------------------------

The basic engine setup is easy to do, and enables the core of
OpenPathSampling's functionality. However, some additional functionality can
be added with very little effort, for engines where this might be relevant. 

Shooting Point Velocity Modification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The original descriptions of TPS used two-way shooting with a ":math:`\Delta
p`" velocity modification. That is, the direction of the velocity on one
particle was changed.

OPS has the ability to do this using the
:class:`.VelocityDirectionModifier`, but this procedure becomes extremely
complicated in condensed phase (periodic) systems with holonomic
constraints. Therefore, OPS only allows this for engines where the snapshots
have the ``n_degrees_of_freedom`` feature. In addition, the engine can tell
OPS not to bother removing linear momentum (e.g., if the engine represents a
particle rolling on a surface, like our toy engine) by setting
``engine.ignore_linear_momentum = True``.

MDTraj support
~~~~~~~~~~~~~~

Several aspects of `MDTraj <http://mdtraj.org>`_ support require an
:class:`mdtraj.Topology` object. The developer can add these by adding a
property/attribute called ``mdtraj_topology`` to their engine. With this,
you automatically get:

* trajectory conversion to :class:`mdtraj.Trajectory`, which enables output
  to various file formats, and integration with other projects such as
  nglview for trajectory visualization.
* ??? TODO: can we make MDTraj CV support automated as well? Check for topol
  on first run? ???



