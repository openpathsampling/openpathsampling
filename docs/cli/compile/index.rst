The Compile Command
===================

The ``compile`` command takes input files and creates OPS object databases
from them. To a large extent, this plays the same role as a traditional input
file for a simulation program: you set the values for certain keywords in a
human-readable file, and the program can create a simulation based on that.

However there are a few significant differences that make the ``compile``
command more flexible and more powerful than a traditional input file:

* **The nested input structure is a better representation of the
  simulation.** This helps provide better context for which keywords are
  relevant, as well as making it more flexible in practical use.
  Additionally, this makes it easier to transition from working with input
  files to writing Python scripts with the library.
* **The input file format allows re-use of named objects within a single
  file.** This not only saves time by avoiding repeated inputs, but because
  the ...
* **Compiled object databases can be re-used.** This is a major difference
  with transition input files. Instead of an input file being specific to a
  single simulation run, the OPS input file creates an object database, and
  that object database can be used for multiple simulation of different
  types. This is also essential for ensuring reproducibility and tracking
  provenance of simulation generated data.


.. toctree::
   :maxdepth: 2

   file_format
   inputs
