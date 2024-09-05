from abc import ABC, abstractmethod
import os
import pathlib
import shutil

import logging
_logger = logging.getLogger(__name__)



class StorageHandler(ABC):
    """Abstract treatment of a key-value-like file/object store.

    This is generally assuming file-based semantics. This may typically mean
    putting things into a temporary directory. This is particularly focused
    on checkpointing, where we will copy the data to put it in a zip file.
    """
    @abstractmethod
    def store(self, storage_label, source_path):
        raise NotImplementedError()

    @abstractmethod
    def load(self, storage_label, target_path):
        raise NotImplementedError()

    @abstractmethod
    def delete(self, storage_label):
        raise NotImplementedError()

    @abstractmethod
    def __contains__(self, storage_label):
        raise NotImplementedError()

    def transfer(self, storage_label, source_path):
        """Transfer a file to the storage label from the source path.

        In some cases, this can be made faster than store followed by
        os.remove, so this method can be overridden. (Example: moving ona
        file system is faster than copying.)
        """
        self.store(storage_label, source_path)
        os.remove(source_path)

    @abstractmethod
    def list_directory(self, storage_label):
        raise NotImplementedError()


class LocalFileStorageHandler(StorageHandler):
    """Concrete implementation of StorageHandler for local files.
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
        os.remove(self.root / storage_label)

    def __contains__(self, storage_label):
        return (self.root / storage_label).exists()

    def transfer(self, storage_label, source_path):
        shutil.move(source_path, self.root / storage_label)

    def list_directory(self, storage_label):
        return [str(p.relative_to(self.root))
                for p in (self.root / storage_label).iterdir()]


class MemoryStorageHandler(StorageHandler):
    """Useful in testing."""
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

        return [key for key in self._data
                if str(key).startswith(str(storage_label))]
