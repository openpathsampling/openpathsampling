from openpathsampling.utils.storage_interfaces import MemoryStorageInterface
from openpathsampling.checkpointing import *

import openpathsampling as paths

import pytest

import pathlib
import zipfile
import tempfile
import shutil

class TestEmptyCheckpointer:
    def setup_method(self):
        self.checkpoint = EmptyCheckpointer()

    def teardown_method(self):
        from openpathsampling.experimental.storage import monkey_patches
        global paths
        paths = monkey_patches.unpatch(paths, with_old_cvs=False)

    def test_context_manager(self):
        with self.checkpoint as cpt:
            assert cpt == self.checkpoint

    def test_load_checkpoint(self):
        assert self.checkpoint.load_checkpoint() == (None, None)

    def test_save_checkpoint(self, tmp_path):
        localfile = tmp_path / "localfile"
        with open(localfile, mode='w') as f:
            f.write("localfile contents")

        self.checkpoint.save_checkpoint({"foo": "bar"}, {"local": localfile})

    def test_next_context(self):
        assert isinstance(self.checkpoint.next_context(), EmptyCheckpointer)

    def test_delete_checkpoint(self):
        self.checkpoint.save_checkpoint({"foo": "bar"})
        self.checkpoint.delete_checkpoint()

    def test_engine_extras(self):
        assert self.checkpoint.engine_extras == []


class CheckpointBreak(Exception):
    """An exception to trigger a checkpoint failure.

    Takes the place of any exception, but we can catch this one in testing.
    """

class TestCheckpointer:
    def setup_method(self):
        storage_handler = MemoryStorageInterface()
        self.checkpoint = Checkpointer(storage_handler)

    def teardown_method(self):
        from openpathsampling.experimental.storage import monkey_patches
        global paths
        paths = monkey_patches.unpatch(paths, with_old_cvs=False)

    def test_setup_teardown(self):
        try:
            assert self.checkpoint._tempdir_manager is None
            tempdir = self.checkpoint._setup_tempdir()
            assert tempdir.exists()
            assert tempdir.is_dir()
            assert self.checkpoint._tempdir_manager is not None
            self.checkpoint._teardown_tempdir(None, None, None)
            assert self.checkpoint._tempdir_manager is None
            assert not tempdir.exists()
        except:
            # kind of naughty to not pass the right stuff to exit
            self.checkpoint_tempdir_manager.__exit__(None, None, None)
            raise

    def test_context_manager_completes(self):
        # with no parent attribute, this will delete the checkpoint at the
        # end and also delete the temporary directory
        store = self.checkpoint.storage_handler
        with self.checkpoint as cpt:
            # TODO: patch test that the temporary directory existed
            assert self.checkpoint._tempdir_manager is None
            cpt.save_checkpoint({"foo": "bar"})
            assert self.checkpoint._tempdir_manager is None
            assert set(store._data) == {pathlib.Path("check.zip")}

        # checkpoint data has been deleted
        assert set(store._data) == set()

    def test_context_manager_errors(self):
        # this will leave the checkpoint behind, but delete the temporary
        # directory
        store = self.checkpoint.storage_handler
        try:
            with self.checkpoint as cpt:
                assert self.checkpoint._tempdir_manager is None
                cpt.save_checkpoint({"foo": "bar"})
                assert self.checkpoint._tempdir_manager is None
                assert set(store._data) == {pathlib.Path("check.zip")}
                raise CheckpointBreak()
        except CheckpointBreak:
            pass

        assert set(store._data) == {pathlib.Path("check.zip")}

    def test_context_manager_cleans_up_children(self):
        store = self.checkpoint.storage_handler
        with self.checkpoint as cpt:
            with cpt.next_context() as cpt1:
                assert self.checkpoint.children == [cpt1]
                cpt1.save_checkpoint({"foo": "bar"})
                assert set(store._data) == {pathlib.Path("0/check.zip")}

            # not cleaned up when we leave the child context
            assert set(store._data) == {pathlib.Path("0/check.zip")}

        # cleaned up when we leave the parent context
        assert set(store._data) == set()

    def test_context_manager_with_parent_keeps_checkpoints(self):
        store = MemoryStorageInterface()
        checkpoint = Checkpointer(store, parent="foo")
        with checkpoint as cpt:
            # TODO: patch test that the temporary directory existed
            assert self.checkpoint._tempdir_manager is None
            cpt.save_checkpoint({"foo": "bar"})
            assert self.checkpoint._tempdir_manager is None
            assert set(store._data) == {pathlib.Path("check.zip")}

        # this time, exit didn't remove checkpoint
        assert set(store._data) == {pathlib.Path("check.zip")}

    def test_sequential_child_contexts(self):
        store = self.checkpoint.storage_handler
        with self.checkpoint as cpt:
            with cpt.next_context() as cpt1:
                cpt1.save_checkpoint({'foo': 'bar'})

            with cpt.next_context() as cpt2:
                cpt2.save_checkpoint({'baz': 'qux'})

            assert set(store._data) == {pathlib.Path("0/check.zip"),
                                        pathlib.Path("1/check.zip")}

        assert set(store._data) == set()

    def test_next_context(self):
        cpt = self.checkpoint.next_context()
        assert cpt.context == pathlib.Path("0")
        assert cpt.storage_handler == self.checkpoint.storage_handler
        assert cpt.frequency == self.checkpoint.frequency
        assert cpt.labeler == self.checkpoint.labeler
        assert cpt.engine_extras == self.checkpoint.engine_extras
        assert cpt.parent == self.checkpoint

    @pytest.mark.parametrize('n_files', [0, 1, 2])
    def test_save_load_checkpoint_round_trip(self, tmp_path, n_files):
        # this is the only unit test that checks with additional files
        input_files = {}
        for suffix in ['a', 'b', 'c'][:n_files]:
            localfile = tmp_path / f"localfile_{suffix}"
            with open(localfile, mode='w') as f:
                f.write(f"localfile_{suffix} contents")

            input_files[f"local_{suffix}"] = localfile

        input_data = {"foo": "bar"}
        with self.checkpoint as cpt:
            with cpt.next_context() as cpt1:
                cpt1.save_checkpoint(input_data, input_files)

            data, files = cpt1.load_checkpoint()
            assert data == input_data
            assert len(files) == n_files
            assert list(files) == list(input_files)
            for key in input_files:
                # files are different
                assert input_files[key] != files[key]

                # contents are the same
                with open(input_files[key], mode='r') as f:
                    input_contents = f.read()
                with open(files[key], mode='r') as f:
                    output_contents = f.read()

                assert input_contents == output_contents

    def test_delete_checkpoint(self):
        store = self.checkpoint.storage_handler
        self.checkpoint.save_checkpoint({"foo": "bar"})
        assert set(store._data) == {pathlib.Path("check.zip")}
        self.checkpoint.delete_checkpoint()
        assert set(self.checkpoint.storage_handler._data) == set()


class CheckpointLifecycleHarness:
    """
    Mix-in for testing entire lifecycle for checkpoints.

    The ``test_checkpoint_lifecycle`` method is essentially an example of
    use (although we do multiple checkpoints in one functions, instead of
    them being in separate functions). If starts with an outer checkpoint,
    which nests inside it two sequential checkpoints. Different subclasses
    can use different setups of the self.store and self.checkpoint in order
    to test different code paths.
    """
    def setup_method(self):
        self.store = self._setup_storage()
        self.checkpoint = Checkpointer(self.store, context="0")

    def teardown_method(self):
        from openpathsampling.experimental.storage import monkey_patches
        global paths
        paths = monkey_patches.unpatch(paths, with_old_cvs=False)

    @staticmethod
    def _setup_storage():
        return MemoryStorageInterface()

    def test_checkpoint_lifecycle(self, tmp_path):
        self.checkpoint_lifecycle(tmp_path)

    def checkpoint_lifecycle(self, tmp_path):
        with self.checkpoint as cpt1:
            data, files = cpt1.load_checkpoint()
            if data is None:
                foo = 1
                cpt1.save_checkpoint({'foo': foo})
                self.assert_after_cpt1_save()
                assert cpt1._tempdir_manager is None
                tmpdir1 = None
            else:
                foo = data['foo']
                self.assert_after_cpt1_load()
                assert cpt1._tempdir_manager is not None
                tmpdir1 = pathlib.Path(cpt1._tempdir_manager.name)

            # after checkpoint 1: state independent of mechanism
            assert not files   # can be None or empty list
            assert foo == 1
            assert pathlib.Path("0/check.zip") in self.store._data

            # run checkpoint 2
            with cpt1.next_context() as cpt2:
                data, files = cpt2.load_checkpoint()
                localdir = tmp_path / "2"
                file2 = localdir / "file2.txt"
                localdir.mkdir()
                if data is None:
                    bar = 2
                    self.make_file_with_contents(file2, "file2 contents")
                    files = {'file2': file2}

                    cpt2.save_checkpoint({'bar': bar}, files)
                    self.assert_after_cpt2_save()
                    assert cpt2._tempdir_manager is None
                    tmpdir2 = None
                else:
                    bar = data['bar']
                    shutil.move(files['file2'], file2)
                    self.assert_after_cpt2_load()
                    assert cpt2._tempdir_manager is not None
                    tmpdir2 = pathlib.Path(cpt2._tempdir_manager.name)

                # after checkpoint 2: state independent of mechanism
                assert bar == 2
                with open(file2, mode='r') as f:
                    assert f.read() == "file2 contents"

                assert pathlib.Path("0/check.zip") in self.store._data
                assert pathlib.Path("0/0/check.zip") in self.store._data

            # after checkpoint 2 exit: independent of mechanism
            if tmpdir2:
                assert not tmpdir2.exists()

            # run checkpoint 3
            with cpt1.next_context() as cpt3:
                data, files = cpt3.load_checkpoint()
                localdir = tmp_path / "3"
                localdir.mkdir()
                file31 = localdir / "file31.txt"
                file32 = localdir / "file32.txt"
                if data is None:
                    baz = 3
                    self.make_file_with_contents(file31, "file31 contents")
                    self.make_file_with_contents(file32, "file32 contents")
                    files = {"file31": file31, "file32": file32}
                    cpt3.save_checkpoint({'baz': baz}, files)
                    self.assert_after_cpt3_save()
                    assert cpt3._tempdir_manager is None
                    tmpdir3 = None
                else:
                    baz = data['baz']
                    shutil.move(files['file31'], file31)
                    shutil.move(files['file32'], file32)
                    self.assert_after_cpt3_load()
                    assert cpt3._tempdir_manager is not None
                    tmpdir3 = pathlib.Path(cpt3._tempdir_manager.name)

                result = self.result_func(foo, bar, baz)
                assert result == 6  # 1 + 2 + 3 (unless we errored)

                # after checkpoint 3: regardless of mechanism
                assert baz == 3
                with open(file31, mode='r') as f:
                    assert f.read() == "file31 contents"
                with open(file32, mode='r') as f:
                    assert f.read() == "file32 contents"

                assert set(self.store._data) == {
                    pathlib.Path("0/check.zip"),
                    pathlib.Path("0/0/check.zip"),
                    pathlib.Path("0/1/check.zip"),
                }

            # after checkpoint 3 exit: independent of mechanism
            if tmpdir3:
                assert not tmpdir3.exists()

        # after checkpoint 1 exit: independent of mechanism
        if tmpdir1:
            assert not tmpdir1.exists()

        # the checkpoints have been cleared after root exits
        assert set(self.store._data) == set()

    def result_func(self, foo, bar, baz):
        return foo + bar + baz

    @staticmethod
    def make_file_with_contents(path, contents):
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, mode='w') as f:
            f.write(contents)

    @staticmethod
    def assert_file_contents(filename, expected):
        with open(filename, mode='r') as f:
            assert f.read() == expected

    @staticmethod
    def assert_zip_directory(zipfilename, expected):
        with zipfile.ZipFile(zipfilename, mode='r') as zipf:
            assert set(zipf.namelist()) == set(expected)

    @staticmethod
    def assert_zip_file_contents(zipfilename, arcname, expected):
        with zipfile.ZipFile(zipfilename, mode='r') as zipf:
            zipcontents = zipf.read(arcname)
            assert zipcontents == expected

    def assert_after_cpt1_load(self):
        raise NotImplementedError("This should not be called")

    def assert_after_cpt1_save(self):
        raise NotImplementedError("This should not be called")

    def assert_after_cpt2_load(self):
        raise NotImplementedError("This should not be called")

    def assert_after_cpt2_save(self):
        raise NotImplementedError("This should not be called")

    def assert_after_cpt3_load(self):
        raise NotImplementedError("This should not be called")

    def assert_after_cpt3_save(self):
        raise NotImplementedError("This should not be called")


class TestStandardLifecycle(CheckpointLifecycleHarness):
    """Standard lifecycle where everything goes well.

    At the end of this, the checkpoint data should be deleted.  All
    checkpoint saves are encounted; all loads are forbidden.
    """
    def assert_after_cpt1_save(self):
        pass

    def assert_after_cpt2_save(self):
        pass

    def assert_after_cpt3_save(self):
        pass


class TestErrorOccursLifecycle(CheckpointLifecycleHarness):
    """Lifecycle where we have a failure while running.

    At the end of this, the checkpoint data should remain. Error occurs
    after all the checkpoints, so we should have checkpoint data from them
    all. Runs all the save methods; all loads should be forbidden.
    """
    def test_checkpoint_lifecycle(self, tmp_path):
        try:
            self.checkpoint_lifecycle(tmp_path)
        except CheckpointBreak:
            assert set(self.store._data) == {
                pathlib.Path("0/check.zip"),
                pathlib.Path("0/0/check.zip"),
                pathlib.Path("0/1/check.zip"),
            }
        else:
            raise AssertionError("Did not cause checkpoint break")

    def result_func(self, foo, bar, baz):
        raise CheckpointBreak()

    def assert_after_cpt1_save(self):
        pass

    def assert_after_cpt2_save(self):
        pass

    def assert_after_cpt3_save(self):
        pass


class TestRecoveryLifecycle(CheckpointLifecycleHarness):
    """Lifecycle where we recover from a checkpoint.

    This starts with checkpoint data in place, and then we delete it all at
    the end. This should follow all the load methods; the saves are
    forbidden.
    """
    @staticmethod
    def _setup_storage():
        store = MemoryStorageInterface()

        with tempfile.TemporaryDirectory() as tmpdir:
            td = pathlib.Path(tmpdir)
            # add checkpoint 1 info
            cpt1 = Checkpointer(store, context="0")
            cpt1.save_checkpoint({"foo": 1})

            # add checkpoint 2 info
            file2 = td / "2/file2.txt"
            file2.parent.mkdir()
            with open(file2, mode='w') as f:
                f.write("file2 contents")
            cpt2 = Checkpointer(store, context="0/0")
            cpt2.save_checkpoint({"bar": 2}, {"file2": file2})

            # add checkpoint 3 info
            file31 = td / "3/file31.txt"
            file31.parent.mkdir()
            with open(file31, mode='w') as f:
                f.write("file31 contents")

            file32 = td / "3/file32.txt"
            with open(file32, mode='w') as f:
                f.write("file32 contents")

            cpt3 = Checkpointer(store, context="0/1")
            cpt3.save_checkpoint({"baz": 3}, {"file31": file31,
                                              "file32": file32})

        assert set(store._data) == {
            pathlib.Path("0/check.zip"),
            pathlib.Path("0/0/check.zip"),
            pathlib.Path("0/1/check.zip"),
        }

        return store

    def assert_after_cpt1_load(self):
        pass

    def assert_after_cpt2_load(self):
        pass

    def assert_after_cpt3_load(self):
        pass
