OpenPathSampling
================

*Read, write and analyze MD trajectories with only a few lines of Python code.*

MDTraj is a python library that allows users to manipulate `molecular dynamics (MD) <http://en.wikipedia.org/wiki/Molecular_dynamics>`_ trajectories.
Extensive :ref:`trajectory analysis <analysis>` routines are implemented. With MDTraj,
you can

 - Read and write from **every MD format imaginable** (``pdb``, ``xtc``, ``trr``,
   ``dcd``, ``binpos``, ``netcdf``, ``mdcrd``, ``prmtop``, ...)
 - Run **blazingly fast** RMSD calculations (4x the speed of the original
   `Theobald QCP <http://theobald.brandeis.edu/qcp/>`_).
 - Use tons of :ref:`analysis <analysis>` functions like bonds/angles/dihedrals,
   hydrogen bonding identification, secondary structure assignment, NMR observables.
 - **Lightweight API**, with a focus on **speed** and vectorized operations.

The library also ships with a flexible command-line application for converting
trajectories between formats. When you install MDTraj, the script will be
installed under the name ``mdconvert``.

.. raw:: html

  <div>
      <h2 style="display: inline; float:left; margin-left:5em">
          <a href="https://github.com/choderalab/opentis/">
          Download the Code</a>
      </h2>
      <h2 style="display: inline; float:right; margin-right:5em">
          <a href="examples/index.html">
          See it in Action</a>
      </h2>
      <div style="clear:both"></div>
      <h2 style="display: inline; float:left; margin-left:7em">
          <a href="http://discourse.openpathsampling.org/">
          Get Help</a>
      </h2>
      <h2 style="display: inline; float:right; margin-right:7em">
          <a href="https://github.com/choderalab/opentis/issues">
          Get Involved</a>
      </h2>
  </div>


.. raw:: html

   <div style="display:none">


--------------------------------------------------------------------------------

Documentation
-------------

.. toctree::
   :maxdepth: 1

   getting_started
   examples/index
   whatsnew
   faq
   Discussion Forums <http://discourse.openpathsampling.org>

API Reference
-------------

.. toctree::
   :maxdepth: 1

   ensemble


Developing
----------

.. toctree::
   :maxdepth: 1

   netcdf_format

--------------------------------------------------------------------------------

.. raw:: html

   </div>

License
-------
OpenPathSampling is licensed under ???