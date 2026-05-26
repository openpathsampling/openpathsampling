from openpathsampling.utils.storage_interfaces import *
import pytest
import tempfile

class StorageInterfaceTest:
    def setup_method(self):
        self.tmpdir_manager = tempfile.TemporaryDirectory()
        self.tmpdir = pathlib.Path(self.tmpdir_manager.__enter__())
        self.localdir = self.tmpdir / "local"
        self.localdir.mkdir()

        # unstored local file
        self.localfile = self.localdir / "localfile"
        with open(self.localfile, mode='w') as f:
            f.write("localfile contents")

        self._initialize_interface()

    def _initialize_interface(self):
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
        self.interface.load("prestored", target_file)
        assert target_file.exists()
        with open(target_file, mode='r') as f:
            assert f.read() == "prestored contents"

    def test_transfer(self):
        raise NotImplementedError()

    def test_transfer_directory(self):
        raise NotImplementedError()

    def test_list_directory(self):
        expected = {"nested/nest_prestored", "nested/deeply/prestored"}
        assert set(self.interface.list_directory("nested")) == expected

    def test_contains(self):
        assert "prestored" in self.interface
        assert "nested/nest_prestored" in self.interface
        assert "nested" not in self.interface
        assert "nonexistent" not in self.interface

    def test_store_existing_no_force(self):
        with pytest.raises(FileExistsError):
            self.interface.store("prestored", self.localfile)

    def test_store_existing_force(self):
        self.interface.store("prestored", self.localfile, force=True)
        assert self.interface.load_bytes("prestored") == b"localfile contents"

    def test_store_bytes(self):
        self.interface.store_bytes("bytes", b"bytes contents")
        assert "bytes" in self.interface
        assert self.interface.load_bytes("bytes") == b"bytes contents"

    def test_store_bytes_existing_no_force(self):
        with pytest.raises(FileExistsError):
            self.interface.store_bytes("prestored", b"new")

    def test_load_bytes(self):
        assert self.interface.load_bytes("prestored") == b"prestored contents"

    def test_open_read(self):
        with self.interface.open("prestored", mode="rb") as f:
            assert f.read() == b"prestored contents"

    def test_open_write(self):
        with self.interface.open("opened", mode="wb") as f:
            f.write(b"opened contents")
        assert self.interface.load_bytes("opened") == b"opened contents"

    def test_open_write_exception_does_not_commit(self):
        with pytest.raises(RuntimeError):
            with self.interface.open("failed-open", mode="wb") as f:
                f.write(b"failed contents")
                raise RuntimeError("fail")
        assert "failed-open" not in self.interface

    def test_open_text_mode_error(self):
        with pytest.raises(ValueError, match="binary"):
            with self.interface.open("prestored", mode="r"):
                pass

    def test_as_path_read(self):
        with self.interface.as_path("prestored", mode="rb") as path:
            with open(path, mode="rb") as f:
                assert f.read() == b"prestored contents"

    def test_as_path_write(self):
        with self.interface.as_path("as-path", mode="wb") as path:
            with open(path, mode="wb") as f:
                f.write(b"as-path contents")
        assert self.interface.load_bytes("as-path") == b"as-path contents"

    def test_as_path_atomic_exception_does_not_commit(self):
        with pytest.raises(RuntimeError):
            with self.interface.as_path("failed-path", mode="wb",
                                        atomic=True) as path:
                with open(path, mode="wb") as f:
                    f.write(b"failed contents")
                raise RuntimeError("fail")
        assert "failed-path" not in self.interface


class TestLocalFileStorageInterface(StorageInterfaceTest):
    def _initialize_interface(self):
        root = self.tmpdir / "stored"
        self.interface = LocalFileStorageInterface(root)
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
        stored_file = self.interface.root / "foo"
        assert not stored_file.exists()
        self.interface.store("foo", self.localfile)
        assert stored_file.exists()
        assert self.localfile.exists()
        with open(stored_file, mode='r') as f:
            assert f.read() == "localfile contents"

    def test_delete(self):
        assert (self.interface.root / 'prestored').exists()
        assert 'prestored' in self.interface
        self.interface.delete('prestored')
        assert "prestored" not in self.interface
        assert not (self.interface.root / 'prestored').exists()

    def test_delete_directory(self):
        nested_file = "nested/nest_prestored"
        assert (self.interface.root / nested_file).exists()
        assert nested_file in self.interface
        with pytest.raises(ValueError, match="is a directory"):
            self.interface.delete('nested')

    def test_transfer(self):
        stored_target = self.interface.root / "foo"
        assert not stored_target.exists()
        self.interface.transfer("foo", self.localfile)
        assert stored_target.exists()
        assert not self.localfile.exists()
        with open(stored_target, mode='r') as f:
            assert f.read() == "localfile contents"

    def test_transfer_nested(self):
        stored_target = self.interface.root / "nested_transfer/foo"
        assert not stored_target.exists()
        self.interface.transfer("nested_transfer/foo", self.localfile)
        assert stored_target.exists()
        assert not self.localfile.exists()
        with open(stored_target, mode='r') as f:
            assert f.read() == "localfile contents"

    def test_as_path_direct_write_uses_final_path(self):
        stored_target = self.interface.root / "direct/foo"
        with self.interface.as_path("direct/foo", mode="wb",
                                    atomic=False) as path:
            assert path == stored_target
            assert not path.exists()
            with open(path, mode="wb") as f:
                f.write(b"direct contents")
        assert stored_target.read_bytes() == b"direct contents"

    def test_as_path_atomic_write_stages_then_commits(self):
        stored_target = self.interface.root / "atomic/foo"
        with self.interface.as_path("atomic/foo", mode="wb",
                                    atomic=True) as path:
            assert path != stored_target
            assert path.parent == stored_target.parent
            with open(path, mode="wb") as f:
                f.write(b"atomic contents")
            assert not stored_target.exists()
        assert stored_target.read_bytes() == b"atomic contents"

    def test_transfer_directory(self):
        source_dir = self.localdir / "directory"
        source_dir.mkdir(exist_ok=True, parents=True)
        subfile = source_dir / "file"
        with open(subfile, mode='w') as f:
            f.write("directory/file contents")

        assert "directory/file" not in self.interface
        assert subfile.exists()

        with pytest.raises(ValueError, match="is a directory"):
            self.interface.transfer("directory", source_dir)

    def test_list_directory_not_directory(self):
        with pytest.raises(ValueError, match="is not a directory"):
            self.interface.list_directory("prestored")


class TestMemoryStorageInterface(StorageInterfaceTest):
    def _initialize_interface(self):
        self.interface = MemoryStorageInterface()
        data = {
            "prestored": b"prestored contents",
            "nested/nest_prestored": b"nested prestored contents",
            "nested/deeply/prestored": b"deeply nested",
        }
        self.interface._data = data

    def test_store(self):
        stored_target = "foo"
        assert stored_target not in self.interface._data
        self.interface.store(stored_target, self.localfile)
        assert stored_target in self.interface._data
        assert self.localfile.exists()
        assert self.interface._data[stored_target] == b"localfile contents"

    def test_transfer(self):
        stored_target = "foo"
        assert stored_target not in self.interface._data
        self.interface.transfer(stored_target, self.localfile)
        assert stored_target in self.interface._data
        assert not self.localfile.exists()
        assert self.interface._data[stored_target] == b"localfile contents"

    def test_delete(self):
        assert "prestored" in self.interface._data
        self.interface.delete("prestored")
        assert "prestored" not in self.interface._data

    def test_delete_directory(self):
        with pytest.raises(KeyError):
            self.interface.delete("nested")

    def test_transfer_directory(self):
        source_dir = self.localdir / "directory"
        source_dir.mkdir(exist_ok=True, parents=True)
        subfile = source_dir / "file"
        with open(subfile, mode='w') as f:
            f.write("directory/file contents")

        assert "directory/file" not in self.interface._data
        assert subfile.exists()

        with pytest.raises(ValueError, match="is a directory"):
            self.interface.transfer("directory", source_dir)

    def test_list_root_directory(self):
        expected = {
            "prestored",
            "nested/nest_prestored",
            "nested/deeply/prestored"
        }
        root = pathlib.Path("")
        assert set(self.interface.list_directory(root)) == expected
