import pathlib
import contextlib
import logging

import numpy as np

_logger = logging.getLogger(__name__)


class TrajectoryWriter:
    """Base class for tools to write trajectories to extenral files.

    This is essentially a wrapper for a function to write the trajectory to
    a file. The function should take two arguments: the trajectory to write,
    and the filename to write to.

    We use an object-oriented approach here so that initialization can
    inlude arbitrary parameters, but the interface used by other code is
    consistent. Additionally, :meth:`__call__` is can handle some standard
    error checking.
    """
    def __call__(self, trajectory, filename, force=False):
        if not force and pathlib.Path(filename).exists():
            raise FileExistsError(f"File {filename} already exists")

        self._write(trajectory, filename)

    @property
    def ext(self):
        """The file extension used by this writer."""
        raise NotImplementedError()

    def _write(self, trajectory, filename):
        """Write the trajectory to the file.

        Parameters
        ----------
        trajectory : openpathsampling.Trajectory
            the trajectory to save
        filename : str
            the name of the file to save to
        """
        raise NotImplementedError()

class SimStoreTrajectoryWriter(TrajectoryWriter):
    """Trajectory writer that uses the OpenPathSampling storage format.

    This is the default trajectory writer, since all engines should be able
    to use it.
    """
    @property
    def ext(self):
        return "db"

    def _write(self, trajectory, filename):
        from openpathsampling.experimental.storage import Storage
        from openpathsampling.experimental.storage.monkey_patches import (
            _IS_PATCHED_SAVING
        )
        if not _IS_PATCHED_SAVING:
            raise RuntimeError("SimStoreTrajectoryWriter requires the "
                               "monkey-patch to be active")
        storage = Storage(filename, mode='w')
        storage.save(trajectory)
