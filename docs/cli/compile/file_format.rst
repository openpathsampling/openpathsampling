Compile input file: Overview
============================

Input files to the ``compile`` command can be in either YAML or JSON format.
Files in YAML format should have the extension ``.yaml`` or ``.yml`` and
files in JSON format should have the extension ``.json`` or ``.jsn``. Most
examples here will use YAML format, although this page will include
equivalent versions in JSON.

The advantages for the YAML format are that it generally requires less
typing, and it allows comments. The main advantage of JSON is that it is
easier to see the structure of mappings and lists than in YAML. Which you
use is personal preference -- the ``compile`` command reads both equally.


Terminology
-----------

There are a few terms that we will use 

* **keyword**: Like most input file formats, an import part of OPS compile
  files is mapping keywords that identify a property to a user-defined
  value for that property.
* **top-level keyword**: Only certain keywords are allowed at the outermost
  level of the nested structure. These keywords are the top-level keywords.
  The top-level keywords in OPS are: ``engines``, ``cvs``, ``volumes``,
  ``states``, ``networks``, and ``schemes``.

Input files for OPS are structures as nested key-value pairs. Each object
that you will store in the database has certain keywords that are needed
define it, and the values you assign to those keywords depend on your
system.


Example: Basics and evaluated expressions
-----------------------------------------

Let's use an OpenMM engine as a simple example. Since this is an engine,
we'll put it under the ``engines`` top-level keyword.  From the
:ref:`documentation on the OpenMM engine <engine--openmm>`, we know that it
requires several parameters:

.. code:: yaml

    engines:
      # note that this is a list of engines, even though we only have one!
      # this is why there is a '-' before the mapping
      - type: openmm  # must be exactly openmm to be identified correctly
        name: my_engine
        topology: topology.pdb       # note that filenames should be relative
        system: system.xml           # to the directory from which you run
        integrator: integrator.xml   # the compile command
        n_steps_per_frame: 1.0 / 0.002  # we allow math here!
        n_frames_max: 10000

We create this within the ``engines`` top-level keyword. An engine can only
be created in a context where the CLI is expecting engines, and the
top-level keyword gives us a way to create that context.

The ``type`` keyword is required for all engines (and, indeed, most
objects), and it must exactly match the identifier string. Note, however,
that comments are allowed.

We give this object the name ``my_engine``. That name is a free choice for
the user, but since it will be used on the command line, it can be useful to
avoid spaces or other characters that require escaping in the shell.

The ``topology``, ``system``, and ``integrator`` parameters take filenames.
These filenames should be relative to the path from which the
``openpathsampling compile`` command will be invoked.

The final two parameters, ``n_steps_per_frame`` and ``n_frames_max``, are
integers. However, for ``n_steps_per_frame``, we provide a string which will
be treated as an :ref:`evaluated expression <expression_eval>`. In
this example, we imagine that perhaps the integrator's timestep is 0.002 ps
and that the user wants to save a frame of the trajectory every 1.0 ps, so
there are 1.0 / 0.002 = 500 steps per frame. Internally, this will be
converted to an integer according to Python's ``int`` builtin. For
``n_frames_max``, we just give the integer value.

The same code could also be implemented in JSON:

.. code:: json

    {"engines": [{
      "type": "openmm",
      "name": "my_engine",
      "topology": "topology.pdb",
      "system": "system.xml",
      "integrator": "integrator.xml",
      "n_steps_per_frame": "1.0 / 0.002",
      "n_frames_max": 10000,
    }]}

In JSON, the ``n_frames_max`` parameter could either be left as an integer
(as it is) or it could be made a string. However, evaluated expressions,
such as ``n_steps_per_frame``, must be strings.

Example: Re-using a named object, optional arguments
----------------------------------------------------

CVs based on MDTraj require a topology, which can be provided either by a
file or by an engine. Engines are always parsed before CVs, so we can use
the engine that we defined above to provide the topology for our CV. Let's
create a couple of MDTraj-based CVs; see our :ref:`documentation on
compiling MDTraj CVs <cvs--mdtraj>` for more details on the keywords. Note
that this assumes it is in the same file as the engine definition above.

.. code:: yaml

    cvs:
      # we're making 2 CVs; each one starts with a '-' so that YAML knows
      # this is a new item in our list
      - type: mdtraj
        name: dist_AB
        topology: my_engine
        func: compute_distances
        kwargs:
          atom_pairs: [[0, 10]]
        # not included: period_min, period_max. They're not required here.

      - type: mdtraj
        name: phi
        topology: my_engine
        func: compute_dihedrals
        kwargs:
          indices: [[4, 6, 8, 14]]
        period_min: -np.pi
        period_max: np.pi
        # Using evaluated expressions for the period

In both CVs, we give a ``type`` and a ``name``, just as we did for the
engine. In this case, we use the engine's ``name`` for the ``topology``
argument. This is how you re-use objects within a file in input files for
the ``compile`` command.
