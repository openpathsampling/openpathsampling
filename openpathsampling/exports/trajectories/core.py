import pathlib
import contextlib
import logging
from openpathsampling.integration_tools import HAS_MDTRAJ

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
    def _write(self, trajectory, filename):
        # TODO: maybe error if not monkey-patched?
        from openpathsampling.experimental.storage import Storage
        storage = Storage(filename, mode='w')
        storage.save(trajectory)
