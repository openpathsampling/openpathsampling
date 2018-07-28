.. _dev-snapshot-features:

Snapshot Features
=================

The :class:`.Snapshot` in OpenPathSampling represents a point in phase
space; that is, in includes all the information about a single time slice of
the simulation. However, what variables define that phase space can vary
from engine to engine. For example, most molecular dynamics engines include
positions, velocities, and box vectors (to define the periodic unit cell).
However, some (such as OpenMM) may require that these come with explicit
units. Perhaps others won't require all of these -- 2D toy models might not
need box vectors, or dynamics based on pure Monte Carlo might not require
velocities. In addition, other engines may have an extended phase space. For
example, there may be approaches that benefit from saving information about
the quantum wavefunction at each time.

Fundamentally, OPS doesn't care about the exact nature of the snapshot.
Whether it has velocities and so forth only matters at a few specific
points; mainly the :class:`.SnapshotModifier` objects or when exporting to
other formats. OPS *does* need to have access to the full set of fields that
define a snapshot, which we refer to as the snapshot features, but it
doesn't usually care what those features represent. Because of this, it is
relatively easy to create engines that use snapshots with arbitrary sets of
features. The implementation of snapshot features in OPS is designed to
allow you to mix and match features.

.. contents::
    :local:

Using Snapshot Features
-----------------------

The basic idea of snapshot features in OPS is that you attach features to a
:class:`.BaseSnapshot` subclass. Note that the features are attached to the
subclass itself, not an instance thereof. Each engine should implement (at
least) one subclass of a snapshot. Typically, the features are grouped into
files, stored in a ``features`` subdirectory within the engine's directory,
which is made into a package by adding an ``__init__.py`` (potentially
empty).

In the simplest way, a sublass of :class:`.BaseSnapshot` just requires using
the ``@openpathsampling.engines.features.base.attach_features`` decorator,
along with a list of modules within the features package. This allows one to
create the entire class using just a list of the features. For example, we
can write

.. code-block:: python

    from openpathsampling.engines.features.base import attach_features
    from openpathsampling.engines import features
    @attach_features([
        features.velocities,
        features.coordinates,
        features.box_vectors,
        features.engine
    ])
    class MySnapshot(BaseSnapshot):
        """Fully functional snapshot"""
        pass

After that simple bit of code, the snapshot will automatically create all
the relevant code for the snapshot This includes:

* the ability to store these snapshots, based on information provided per
  feature
* information about how to copy the snapshot, including whether each feature
  should be deep copied or shallow copied
* all other properties related to the feature (for example, the default
  ``coordinates`` feature may carry units with it -- the unitless version is
  available in a numpy array accessible as ``snapshot.xyz``. This is
  implemented as a property in the ``coordinates`` module, and is
  automatically included when the ``coordinates`` feature module is attached
  to the snapshot.

One important point: **all snapshots should include the engine feature.**
Several parts of OPS storage and analysis assume that the ``engine``
property is available for any :class:`.Snapshot`.

For many molecular dynamics purposes, the built-in snapshot features will be
sufficient, and you can just attach them to your engine. See the list of
common snapshot features, along with their implementation details, for more
information on the built-in snapshot features. The following subsections
will describe how to create custom snapshot feature modules.

Writing Snapshot Features
-------------------------

Snapshot features are collected into modules that group related features
together. A snapshot feature module may include the following information:

* docstring: A partial numpydoc-style docstring to document the features in
  this module
* feature-defining lists (``variables``, ``minus``, ``numpy``, ``lazy``,
  ``dimensions``): Lists of strings for stored features. See details in the
  section on stored snapshot features.
* ``netcdfplus_init`` method: Method used in initialize storage information
  in the netcdfplus format.
* ``@property`` features: Snapshot features that can be implemented as
  properties. See section on creating property snapshot features for more
  information.

The following subsections will describe how to implement various kinds of
snapshot features within a module. See the list of standard snapshot
features for links to code showing more example implementations.

Creating Non-Stored (Property) Snapshot Features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may want to add properties to snapshots that are not explicitly stored.
For example, the masses of each particle is commonly a constant for a given
instance of an MD engine. However, you might want to access them as
``snapshot.masses``.

To do this with a snapshot feature, you can use the ``@property`` decorator,
just as you would with a property of a class instance. For example, this
might be implemented as

.. code-block:: python

    @property
    def masses(snapshot):
        return snapshot.engine.get_masses()

where we assume that ``engine.get_masses()`` returns the masses. By putting
this in a module called ``masses.py`` and attaching that module as a
feature, the snapshot will automatically have the ``masses``  property.

One common pattern is for features that are the same for all snapshots
created by the engine (such as masses, or number of degrees of freedom) to
use a property of the engine object that is only computed/stored once. The
actual numbers are resident in memory as part of the engine, but are
accessible from any snapshot created by that engine.

Creating Stored Snapshot Features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Features that contain information that should be stored are a bit more
complicated. First, such objects should be registered as "variables" by
including their names in the list of strings in ``variables``.

.. TODO: finish this documentation

Creating Proxy (Container) Snapshot Features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In many cases, we don't want to fully load the information in a snapshot,
such as the coordinates or the velocities. For example, when calculating
something like a histogram of path lengths for a given ensemble, we don't
actually need the coordinates. In order to load snapshots containing
information that is stored, but without loading that information, we use an
extra layer of abstraction called a "container" feature.

One example is the ``statics`` container. This includes the positions and
box vectors for a given snapshot. However, the snapshot can load a pointer
to the statics container without actually loading the positions, thus saving
time and memory.

.. TODO: add details on creating such containers

A similar idea is used for external snapshots, where all data is stored in
an external file. For the implementation of external snapshots, see
documentation on the indirect engine API (coming in version 1.1).

Standard Snapshot Features
--------------------------

In order to help simulation and analysis code to be useful for many engines,
we have some recommended names for snapshot features. By using these names
with the snapshots from your engines, you can automatically gain additional
functionality from other parts of OPS. For example, this enables us to use
the same API when dealing with coordinates whether they are directly
attributes of the snapshot, as with the toy engine, or whether they are
within an additional abstraction layer in a ``statics`` object, as in the
OpenMM engine.

+---------------------+----------------------------------------------------+
| Name                |  Description and implementation examples           |
+=====================+====================================================+
| ``engine``          | :class:`.DynamicsEngine` instance that created     |
|                     | this snapshot. As stored feature in                |
|                     | :mod:`.engines.features.engine`.                   |
+---------------------+----------------------------------------------------+
| ``coordinates``     | Particle positions. Unitted. As stored feature in  |
|                     | :mod:`.engines.features.coordinates`. As property  |
|                     | feature in :mod:`.engines.features.statics`.       |
+---------------------+----------------------------------------------------+
| ``velocities``      | Particle velocities. Unitted. As stored feature    |
|                     | in :mod:`.engines.features.velocities`. As         |
|                     | property feature in                                |
|                     | :mod:`.engines.features.kinetics`.                 |
+---------------------+----------------------------------------------------+
| ``box_vectors``     | Unit cell vectors for the periodic box. Unitted.   |
|                     | As stored feature in                               |
|                     | :mod:`.engines.features.box_vectors` As property   |
|                     | feature in :mod:`.engines.features.statics`.       |
+---------------------+----------------------------------------------------+
| ``statics``         | Combination of coordinates and box vectors. As     |
|                     | container feature in                               |
|                     | :mod:`.engines.features.statics` and               |
|                     | :mod:`.engines.features.shared`.                   |
+---------------------+----------------------------------------------------+
| ``kinetics``        | Container for velocities. As container feature in  |
|                     | :mod:`.engines.features.kinetics` and              |
|                     | :mod:`.engines.features.shared`.                   |
+---------------------+----------------------------------------------------+
| ``xyz``             | Particle positions, without units. As property     |
|                     | feature in :mod:`.engines.features.coordinates`,   |
|                     | :mod:`.engines.features.statics`.                  |
+---------------------+----------------------------------------------------+
| ``masses``          | Particle masses (in actual mass units, not mass    |
|                     | per mole, as used in some engines). Unitted.       |
|                     | As property feature in                             |
|                     | :mod:`.engines.openmm.features.masses`.            |
+---------------------+----------------------------------------------------+
| ``mass_per_mole``   | Particle mass per mole. Used in as mass in some    |
|                     | engines to provide energies in per-mole units.     |
|                     | Unitted. As property feature in                    |
|                     | :mod:`.engines.openmm.features.masses`.            |
+---------------------+----------------------------------------------------+
| |ndof|              | Number of degrees of freedom. Should account for   |
|                     | any constraints (including, e.g., total linear     |
|                     | momentum.) As property feature in                  |
|                     | |openmm_inst_temp|, |toy_inst_temp|.               |
+---------------------+----------------------------------------------------+
| |inst_temp|         | Instantaneous temperature, i.e., the effective     |
|                     | temperature at a point in time based on the        |
|                     | kinetic energy. As property feature in             |
|                     | |openmm_inst_temp|, |toy_inst_temp|.               |
+---------------------+----------------------------------------------------+

.. |ndof| replace:: ``n_degrees_of_freedom``
.. |inst_temp| replace:: ``instantaneous_temperature``
.. |openmm_inst_temp| replace:: :mod:`.engines.openmm.features.instantaneous_temperature`
.. |toy_inst_temp| replace:: :mod:`.engines.toy.features.instantaneous_temperature`
