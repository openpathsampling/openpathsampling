from abc import ABC, abstractmethod
import builtins
import contextlib
import io
import os
import pathlib
import shutil
import tempfile

import logging
_logger = logging.getLogger(__name__)



class StorageInterface(ABC):
    """Abstract treatment of a key-value-like file/object store.

    This is generally assuming file-based semantics. This may typically mean
    putting things into a temporary directory. This is particularly focused
    on checkpointing, where we will copy the data to put it in a zip file.

    The default implementations are intended for single-process export
    workflows. ``force=False`` performs a best-effort existence check, but
    this interface does not provide cross-process locking or distributed
    consistency guarantees.
    """
    @abstractmethod
    def store(self, storage_label, source_path, *, force=False):
        """Store the data in ``source_path`` at key ``storage_label``

        Parameters
        ----------
        storage_label : str
            The key to store the data at.
        source_path : os.PathLike
            The path to the data to store.
        force : bool
            If False, raise FileExistsError when ``storage_label`` exists.
        """
        raise NotImplementedError()

    @abstractmethod
    def load(self, storage_label, target_path):
        """Load the data from ``storage_label`` into file at ``target_path``

        Parameters
        ----------
        storage_label : str
            The key to load the data from.
        target_path : os.PathLike
            The path to store the data at.
        """
        raise NotImplementedError()

    @abstractmethod
    def delete(self, storage_label):
        """Delete key ``storage_label`` from the object store.
        """
        raise NotImplementedError()

    @abstractmethod
    def __contains__(self, storage_label):
        raise NotImplementedError()

    def transfer(self, storage_label, source_path, *, force=False):
        """Transfer a file to the storage label from the source path.

        In some cases, this can be made faster than store followed by
        os.remove, so this method can be overridden. (Example: moving on a
        file system is faster than copying.)
        """
        if pathlib.Path(source_path).is_dir():
            raise ValueError(f"'{source_path}' is a directory, and can't "
                             "be transferred.")
        self.store(storage_label, source_path, force=force)
        os.remove(source_path)

    def store_bytes(self, storage_label, data, *, force=False):
        """Store bytes at key ``storage_label``."""
        fd, tmp_path = tempfile.mkstemp()
        tmp_path = pathlib.Path(tmp_path)
        try:
            with os.fdopen(fd, mode='wb') as f:
                f.write(data)
            self.store(storage_label, tmp_path, force=force)
        finally:
            if tmp_path.exists():
                os.remove(tmp_path)

    def load_bytes(self, storage_label):
        """Load bytes stored at key ``storage_label``."""
        fd, tmp_path = tempfile.mkstemp()
        os.close(fd)
        tmp_path = pathlib.Path(tmp_path)
        try:
            self.load(storage_label, tmp_path)
            with builtins.open(tmp_path, mode='rb') as f:
                return f.read()
        finally:
            if tmp_path.exists():
                os.remove(tmp_path)

    @contextlib.contextmanager
    def open(self, storage_label, mode='rb', *, force=False):
        """Open an in-memory binary file object for ``storage_label``.

        Only binary modes are supported. Write modes commit on successful
        context exit; exceptions leave the storage label unchanged.
        """
        if 'b' not in mode:
            raise ValueError("StorageInterface.open only supports binary "
                             "modes")
        if _mode_writes(mode):
            if not force and storage_label in self:
                raise FileExistsError(f"Storage label {storage_label} "
                                      "already exists")
            fileobj = io.BytesIO()
            yield fileobj
            self.store_bytes(storage_label, fileobj.getvalue(), force=force)
        else:
            fileobj = io.BytesIO(self.load_bytes(storage_label))
            yield fileobj

    @contextlib.contextmanager
    def as_path(self, storage_label, mode='rb', *, suffix=None, force=False,
                atomic=False):
        """Yield a local filesystem path for ``storage_label``.

        This supports libraries that only write to filenames. For non-local
        storage and for local atomic writes, implementations may stage through
        a temporary path and commit on successful context exit.
        """
        suffix = suffix or pathlib.Path(storage_label).suffix
        fd, tmp_path = tempfile.mkstemp(suffix=suffix)
        os.close(fd)
        tmp_path = pathlib.Path(tmp_path)
        write_mode = _mode_writes(mode)
        try:
            if write_mode:
                if not force and storage_label in self:
                    raise FileExistsError(f"Storage label {storage_label} "
                                          "already exists")
            else:
                self.load(storage_label, tmp_path)

            yield tmp_path

            if write_mode:
                self.transfer(storage_label, tmp_path, force=force)
        finally:
            if tmp_path.exists():
                os.remove(tmp_path)

    @abstractmethod
    def list_directory(self, storage_label):
        """List all objects in subdirectories of the given storage label.
        """
        raise NotImplementedError()


def _mode_writes(mode):
    """Return True when a file mode can modify storage contents."""
    return any(flag in mode for flag in ('w', 'a', 'x', '+'))


class LocalFileStorageInterface(StorageInterface):
    """Concrete implementation of StorageInterface for local files.

    Parameters
    ----------
    root : os.PathLike
        The root directory for the storage interface.
    """
    def __init__(self, root):
        self.root = pathlib.Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.root = self.root.resolve()

    def _local_path(self, storage_label):
        storage_label = pathlib.Path(storage_label)
        if storage_label.anchor:
            raise ValueError("Storage labels must be unanchored relative paths")

        local_path = (self.root / storage_label).resolve(strict=False)
        try:
            local_path.relative_to(self.root)
        except ValueError as exc:
            raise ValueError("Storage labels must resolve inside the storage "
                             "root") from exc
        return local_path

    def _check_can_write(self, storage_label, force):
        local_path = self._local_path(storage_label)
        if local_path.exists():
            if local_path.is_dir():
                raise ValueError(f"Storage label {storage_label} refers to an "
                                 "existing directory")
            if not force:
                raise FileExistsError(f"Storage label {storage_label} "
                                      "already exists")
        return local_path

    def store(self, storage_label, source_path, *, force=False):
        local_path = self._check_can_write(storage_label, force)
        local_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source_path, local_path)

    def load(self, storage_label, target_path):
        target_path = pathlib.Path(target_path)
        target_path.parent.mkdir(parents=True, exist_ok=True)
        local_path = self._local_path(storage_label)
        _logger.debug("Copying file from {str(local_path)} "
                      f"to {str(target_path)}")
        shutil.copyfile(local_path, target_path)

    def delete(self, storage_label):
        obj = self._local_path(storage_label)
        if obj.is_dir():
            raise ValueError(f"'{obj}' is a directory, and can't be "
                             "deleted.")
        else:
            _logger.debug("Deleting file {str(obj)}")
            os.remove(obj)

    def __contains__(self, storage_label):
        expected = self._local_path(storage_label)
        return expected.exists() and not expected.is_dir()

    def transfer(self, storage_label, source_path, *, force=False):
        if pathlib.Path(source_path).is_dir():
            raise ValueError(f"'{source_path}' is a directory, and can't "
                             "be transferred.")
        target_path = self._check_can_write(storage_label, force)
        target_path.parent.mkdir(parents=True, exist_ok=True)
        if force and target_path.exists():
            os.remove(target_path)
        shutil.move(source_path, target_path)

    def store_bytes(self, storage_label, data, *, force=False):
        local_path = self._check_can_write(storage_label, force)
        local_path.parent.mkdir(parents=True, exist_ok=True)
        with builtins.open(local_path, mode='wb') as f:
            f.write(data)

    def load_bytes(self, storage_label):
        with builtins.open(self._local_path(storage_label), mode='rb') as f:
            return f.read()

    @contextlib.contextmanager
    def open(self, storage_label, mode='rb', *, force=False):
        if 'b' not in mode:
            raise ValueError("StorageInterface.open only supports binary "
                             "modes")
        if _mode_writes(mode):
            with super().open(storage_label, mode=mode, force=force) as f:
                yield f
        else:
            local_path = self._local_path(storage_label)
            with builtins.open(local_path, mode=mode) as f:
                yield f

    @contextlib.contextmanager
    def as_path(self, storage_label, mode='rb', *, suffix=None, force=False,
                atomic=False):
        """Yield a local path for ``storage_label``.

        With ``atomic=False`` and a write mode, this yields the final path
        directly for performance and disk-space efficiency. Failed or
        interrupted writes may therefore leave partial contents. Use
        ``atomic=True`` to stage through a temporary file and commit only on
        successful context exit. Neither mode provides cross-process locking.
        """
        local_path = self._local_path(storage_label)
        write_mode = _mode_writes(mode)
        if not atomic:
            if write_mode:
                local_path = self._check_can_write(storage_label, force)
                local_path.parent.mkdir(parents=True, exist_ok=True)
            yield local_path
            return

        if write_mode:
            local_path = self._check_can_write(storage_label, force)
            local_path.parent.mkdir(parents=True, exist_ok=True)
            suffix = suffix or local_path.suffix
            fd, tmp_path = tempfile.mkstemp(
                suffix=suffix,
                dir=str(local_path.parent)
            )
            os.close(fd)
            tmp_path = pathlib.Path(tmp_path)
            try:
                yield tmp_path
                if not force and storage_label in self:
                    raise FileExistsError(f"Storage label {storage_label} "
                                          "already exists")
                os.replace(tmp_path, local_path)
            finally:
                if tmp_path.exists():
                    os.remove(tmp_path)
        else:
            yield local_path

    def list_directory(self, storage_label):
        path = self._local_path(storage_label)
        if not path.is_dir():
            raise ValueError(f"'{path}' is not a directory.")
        return [
            str((pathlib.Path(p[0]) / subp).relative_to(self.root))
            for p in os.walk(path)
            for subp in p[2]
        ]


class MemoryStorageInterface(StorageInterface):
    """In-memory storage interface.

    Useful in testing.
    """
    def __init__(self):
        self._data = {}

    def _check_can_write(self, storage_label, force):
        if not force and storage_label in self:
            raise FileExistsError(f"Storage label {storage_label} "
                                  "already exists")

    def store(self, storage_label, source_path, *, force=False):
        self._check_can_write(storage_label, force)
        with builtins.open(source_path, mode='rb') as f:
            self._data[storage_label] = f.read()

    def load(self, storage_label, target_path):
        target_path = pathlib.Path(target_path)
        target_path.parent.mkdir(parents=True, exist_ok=True)
        with builtins.open(target_path, mode='wb') as f:
            f.write(self._data[storage_label])
        return self._data[storage_label]

    def store_bytes(self, storage_label, data, *, force=False):
        self._check_can_write(storage_label, force)
        self._data[storage_label] = data

    def load_bytes(self, storage_label):
        return self._data[storage_label]

    def delete(self, storage_label):
        del self._data[storage_label]

    def __contains__(self, storage_label):
        return storage_label in self._data

    def list_directory(self, storage_label):
        # special case because the empty path becomes '.' as a string
        if storage_label == pathlib.Path(""):
            storage_label = ""
        elif not storage_label.endswith("/"):
            storage_label += "/"

        return [key for key in self._data
                if str(key).startswith(str(storage_label))]
