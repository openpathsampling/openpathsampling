from openpathsampling.utils.storage_handlers import *
import pytest
import tempfile

class StorageHandlerTest:
    def setup_method(self):
        self.tmpdir_manager = tempfile.TemporaryDirectory()
        self.tmpdir = pathlib.Path(self.tmpdir_manager.__enter__())
        self.localdir = self.tmpdir / "local"
        self.localdir.mkdir()

        # unstored local file
        self.localfile = self.localdir / "localfile"
        with open(self.localfile, mode='w') as f:
            f.write("localfile contents")

        self._initialize_handler()

    def _initialize_handler(self):
        """Subclasses must implement this initialization routine.

        This method must create the following stored objects:

        * key ``prestored``, contents ``prestored contents``
        * key ``nested/nest_prestored``, contents ``nested prestored
          contents``
        * key ``nested/deeply/prestored``, contents ``deeply nested``
        """
        raise NotImplementedError()

    def teardown_method(self):
        self.tmpdir_manager.__exit__(None, None, None)

    def test_store(self):
        raise NotImplementedError()

    def test_delete(self):
        raise NotImplementedError()

    def test_delete_directory(self):
        raise NotImplementedError()

    def test_load(self):
        target_file = self.localdir / "foo"
        assert not target_file.exists()
        self.handler.load("prestored", target_file)
        assert target_file.exists()
        with open(target_file, mode='r') as f:
            assert f.read() == "prestored contents"

    def test_transfer(self):
        raise NotImplementedError()

    def test_transfer_directory(self):
        raise NotImplementedError()

    def test_list_directory(self):
        expected = {"nested/nest_prestored", "nested/deeply/prestored"}
        assert set(self.handler.list_directory("nested")) == expected

class TestLocalFileStorageHandler(StorageHandlerTest):
    def _initialize_handler(self):
        root = self.tmpdir / "stored"
        self.handler = LocalFileStorageHandler(root)
        # pre-stored file
        with open(root / "prestored", mode='w') as f:
            f.write("prestored contents")

        # pre-stored nested files (for directory testing)
        nested_prestored = root / "nested/nest_prestored"
        nested_prestored.parent.mkdir()
        with open(nested_prestored, mode='w') as f:
            f.write("nested prestored contents")

        # deeply nested file
        deeply_prestored = root / "nested/deeply/prestored"
        deeply_prestored.parent.mkdir(parents=True)
        with open(deeply_prestored, mode='w') as f:
            f.write("deeply nested")

    def test_store(self):
        stored_file = self.handler.root / "foo"
        assert not stored_file.exists()
        self.handler.store("foo", self.localfile)
        assert stored_file.exists()
        assert self.localfile.exists()
        with open(stored_file, mode='r') as f:
            assert f.read() == "localfile contents"

    def test_delete(self):
        assert (self.handler.root / 'prestored').exists()
        assert 'prestored' in self.handler
        self.handler.delete('prestored')
        assert "prestored" not in self.handler
        assert not (self.handler.root / 'prestored').exists()

    def test_delete_directory(self):
        nested_file = "nested/nest_prestored"
        assert (self.handler.root / nested_file).exists()
        assert nested_file in self.handler
        with pytest.raises(ValueError, match="is a directory"):
            self.handler.delete('nested')

    def test_transfer(self):
        stored_target = self.handler.root / "foo"
        assert not stored_target.exists()
        self.handler.transfer("foo", self.localfile)
        assert stored_target.exists()
        assert not self.localfile.exists()
        with open(stored_target, mode='r') as f:
            assert f.read() == "localfile contents"

    def test_transfer_directory(self):
        source_dir = self.localdir / "directory"
        source_dir.mkdir(exist_ok=True, parents=True)
        subfile = source_dir / "file"
        with open(subfile, mode='w') as f:
            f.write("directory/file contents")

        assert "directory/file" not in self.handler
        assert subfile.exists()

        with pytest.raises(ValueError, match="is a directory"):
            self.handler.transfer("directory", source_dir)


class TestMemoryStorageHandler(StorageHandlerTest):
    def _initialize_handler(self):
        self.handler = MemoryStorageHandler()
        data = {
            "prestored": b"prestored contents",
            "nested/nest_prestored": b"nested prestored contents",
            "nested/deeply/prestored": b"deeply nested",
        }
        self.handler._data = data

    def test_store(self):
        stored_target = "foo"
        assert stored_target not in self.handler._data
        self.handler.store(stored_target, self.localfile)
        assert stored_target in self.handler._data
        assert self.localfile.exists()
        assert self.handler._data[stored_target] == b"localfile contents"

    def test_transfer(self):
        stored_target = "foo"
        assert stored_target not in self.handler._data
        self.handler.transfer(stored_target, self.localfile)
        assert stored_target in self.handler._data
        assert not self.localfile.exists()
        assert self.handler._data[stored_target] == b"localfile contents"

    def test_delete(self):
        assert "prestored" in self.handler._data
        self.handler.delete("prestored")
        assert "prestored" not in self.handler._data

    def test_delete_directory(self):
        with pytest.raises(KeyError):
            self.handler.delete("nested")

    def test_transfer_directory(self):
        source_dir = self.localdir / "directory"
        source_dir.mkdir(exist_ok=True, parents=True)
        subfile = source_dir / "file"
        with open(subfile, mode='w') as f:
            f.write("directory/file contents")

        assert "directory/file" not in self.handler._data
        assert subfile.exists()

        with pytest.raises(ValueError, match="is a directory"):
            self.handler.transfer("directory", source_dir)
