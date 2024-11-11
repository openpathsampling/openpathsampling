from abc import ABC, abstractmethod
import os
import pathlib
import shutil

import logging
_logger = logging.getLogger(__name__)



class StorageInterface(ABC):
    """Abstract treatment of a key-value-like file/object store.

    This is generally assuming file-based semantics. This may typically mean
    putting things into a temporary directory. This is particularly focused
    on checkpointing, where we will copy the data to put it in a zip file.
    """
    @abstractmethod
    def store(self, storage_label, source_path):
        """Store the data in ``source_path`` at key ``storage_label``

        Parameters
        ----------
        storage_label : str
            The key to store the data at.
        source_path : os.PathLike
            The path to the data to store.
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

    def transfer(self, storage_label, source_path):
        """Transfer a file to the storage label from the source path.

        In some cases, this can be made faster than store followed by
        os.remove, so this method can be overridden. (Example: moving on a
        file system is faster than copying.)
        """
        if pathlib.Path(source_path).is_dir():
            raise ValueError(f"'{source_path}' is a directory, and can't "
                             "be transferred.")
        self.store(storage_label, source_path)
        os.remove(source_path)

    @abstractmethod
    def list_directory(self, storage_label):
        """List all objects in subdirectories of the given storage label.
        """
        raise NotImplementedError()


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

    def store(self, storage_label, source_path):
        local_path = self.root / storage_label
        local_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source_path, self.root / storage_label)

    def load(self, storage_label, target_path):
        target_path.parent.mkdir(parents=True, exist_ok=True)
        local_path = self.root / storage_label
        _logger.debug("Copying file from {str(local_path)} "
                      f"to {str(target_path)}")
        shutil.copyfile(local_path, target_path)

    def delete(self, storage_label):
        obj = self.root / storage_label
        if obj.is_dir():
            raise ValueError(f"'{obj}' is a directory, and can't be "
                             "deleted.")
        else:
            _logger.debug("Deleting file {str(obj)}")
            os.remove(obj)

    def __contains__(self, storage_label):
        expected = self.root / storage_label
        return expected.exists() and not expected.is_dir()

    def transfer(self, storage_label, source_path):
        if pathlib.Path(source_path).is_dir():
            raise ValueError(f"'{source_path}' is a directory, and can't "
                             "be transferred.")
        shutil.move(source_path, self.root / storage_label)

    def list_directory(self, storage_label):
        path = self.root / storage_label
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

    def store(self, storage_label, source_path):
        with open(source_path, mode='rb') as f:
            self._data[storage_label] = f.read()

    def load(self, storage_label, target_path):
        with open(target_path, mode='wb') as f:
            f.write(self._data[storage_label])
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
