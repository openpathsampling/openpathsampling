Trajectory Export from Engines
==============================

Users frequently want trajectories to be exported from Python into format
that they can analyze with other tools. Since the expected trajectory format
can differ between engines, a developer can set a default trajectory writer
using the :meth:`._default_trajectory_writer` method. This method should
return a callable that takes a trajectory, a filename, and the boolean flag
``force``, where the ``force`` flag forces overwrite of an existing file (a
``FileExistsError`` should be raised if the file exists and ``force`` is
``False``. Note that the default writer should be lossless: output from this
should be sufficient to perform an MD restart for, e.g., a shooting move.

We recommend creating this callable as a subclass for the
:class:`.TrajectoryWriter` class. To implement this, you only need to
implement the :meth:`.TrajectoryWriter._write` method, which takes the
trajectory and the filename to write to as arguments. Additional information
that can't be extracted from the trajectory can be passed in the class
``__init__`` method.

If you do not implement :meth:`._default_trajectory_writer`, you will use
the default from the base class, which creates a SimStore file for the
trajectory. OPS should be able to read/write to SimStore for snapshots from
any engine, so this will work, although it may not be convenient to users.
