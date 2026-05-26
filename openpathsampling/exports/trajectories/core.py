import pathlib
import logging

from openpathsampling.utils.storage_interfaces import LocalFileStorageInterface

_logger = logging.getLogger(__name__)


class TrajectoryWriter:
    """Base class for tools to write trajectories to external storage.

    This is essentially a wrapper for a function to write the trajectory to
    a path. Most trajectory libraries require a local filename, so the
    primary storage-aware API uses :meth:`StorageInterface.as_path` to
    provide one.

    We use an object-oriented approach here so that initialization can
    include arbitrary parameters, but the interface used by other code is
    consistent. Additionally, :meth:`write` handles standard error checking
    through the storage interface.
    """
    def __call__(self, trajectory, filename, force=False):
        self.write_file(trajectory, filename, force=force)

    def write(self, trajectory, storage, storage_label, *, force=False,
              atomic=False):
        """Write ``trajectory`` to ``storage_label`` in ``storage``.

        Parameters
        ----------
        trajectory : openpathsampling.Trajectory
            the trajectory to save
        storage : openpathsampling.utils.storage_interfaces.StorageInterface
            destination storage interface
        storage_label : str
            storage key for the trajectory
        force : bool
            overwrite an existing storage label if True
        atomic : bool
            if True, request staged/atomic path handling from storage
        """
        suffix = self._suffix_for_label(storage_label)
        with storage.as_path(storage_label, mode='wb', suffix=suffix,
                             force=force, atomic=atomic) as path:
            self._write_path(trajectory, path)

    def write_file(self, trajectory, filename, *, force=False):
        """Write ``trajectory`` directly to a filesystem path."""
        filename = pathlib.Path(filename)
        storage = LocalFileStorageInterface(filename.parent)
        self.write(trajectory, storage, filename.name, force=force,
                   atomic=False)

    def _suffix_for_label(self, storage_label):
        suffix = pathlib.Path(storage_label).suffix
        if suffix:
            return suffix
        try:
            ext = self.ext
        except NotImplementedError:
            return None
        return f".{ext}" if ext else None

    @property
    def ext(self):
        """The file extension used by this writer."""
        raise NotImplementedError()

    def _write_path(self, trajectory, filename):
        """Write the trajectory to the local filesystem path ``filename``."""
        self._write(trajectory, filename)

    def _write(self, trajectory, filename):
        """Write the trajectory to the local filesystem path ``filename``.

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
