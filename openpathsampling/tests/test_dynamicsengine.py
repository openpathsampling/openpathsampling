from __future__ import absolute_import
from builtins import object
import pytest
import openpathsampling as paths

from .test_helpers import make_1d_traj, raises_with_message_like

try:
    from unittest import mock
except ImportError:
    import mock


class StupidEngine(paths.engines.DynamicsEngine):
    _default_options = {'random_option': False}

    @property
    def bad_property(self):
        obj = object()
        return obj.b

    @property
    def property_recovers(self):
        if not hasattr(self, 'attempted'):
            self.attempted = True
            raise AttributeError("Internal error")
        return self.attempted

    def generate_next_frame(self):
        return make_1d_traj([0.0])[0]

    @property
    def current_snapshot(self):
        return make_1d_traj([0.0])[0]

    @current_snapshot.setter
    def current_snapshot(self, snap):
        pass


class TestDynamicsEngine(object):
    def setup_method(self):
        options = {'n_frames_max' : 100, 'random_option' : True}

        snapshot_dimensions = {
            'n_atoms': 1,
            'n_spatial': 1
        }

        descriptor = paths.engines.snapshot.SnapshotDescriptor.construct(
            snapshot_class=paths.engines.toy.Snapshot,
            snapshot_dimensions=snapshot_dimensions
        )
        self.descriptor = descriptor
        self.engine = paths.engines.DynamicsEngine(options, descriptor)
        self.stupid = StupidEngine(options, descriptor)

    def test_getattr_from_options(self):
        assert self.stupid.random_option is True

    def test_getattr_property_fails(self):
        # py2: "'newobject' object has no attribute 'b'"
        # py3: "'object' object has no attribute 'b'"
        error_string_end = "object' object has no attribute 'b'"
        try:
            self.stupid.bad_property
        except AttributeError as e:
            assert str(e).endswith(error_string_end)

        else:
            raise AssertionError("Did not raise appropriate AttributeError")


    @raises_with_message_like(AttributeError,
                              "Unknown problem occurred in property")
    def test_getattr_property_recovers(self):
        self.stupid.property_recovers

    @raises_with_message_like(AttributeError,
                              "'StupidEngine' has no attribute 'foo'" +
                              ", nor does its options dictionary")
    def test_getattr_does_not_exist(self):
        self.stupid.foo

    def test_getattr_dimension(self):
        assert self.engine.n_atoms == 1
        assert self.engine.n_spatial == 1
        assert self.stupid.n_atoms == 1
        assert self.stupid.n_spatial == 1

    def test_unicode_options(self):
        # regression test: this should NOT raise an error
        options = {'on_nan': u'fail'}
        engine = paths.engines.DynamicsEngine(options, self.descriptor)

    @pytest.mark.parametrize('stoppable', ['first', 'last'])
    def test_generate_multiple_running_conditions(self, stoppable):
        # test to ensure that, when multiple continue conditions are
        # defined, we stop if any of them tell us to stop
        mock_true_false_true = mock.Mock(side_effect=[True, False, True])
        mock_true = mock.Mock(return_value=True)

        init_snap = make_1d_traj([0.0])[0]
        conditions = {'first': [mock_true_false_true, mock_true],
                      'last': [mock_true, mock_true_false_true]}[stoppable]
        # this raises an StopIteration (and then RuntimeError) when it
        # doesn't work
        traj = self.stupid.generate(init_snap, conditions)
        assert len(traj) == 2

    def test_export_trajectory(self, tmp_path):
        global paths
        outfile = tmp_path / "test_traj.db"
        init_snap = make_1d_traj([0.0])[0]
        continue_condition = paths.LengthEnsemble(5).can_append
        traj = self.stupid.generate(init_snap, [continue_condition])
        assert len(traj) == 5
        assert not outfile.exists()
        from openpathsampling.experimental.storage import (
            monkey_patches, Storage
        )
        paths = monkey_patches.monkey_patch_all(paths)
        try:
            self.stupid.export_trajectory(traj, outfile)
            assert outfile.exists()
            storage = Storage(outfile, mode='r')
            assert len(storage.trajectories) == 1
            reloaded = storage.trajectories[0]
            assert len(reloaded) == 5
            assert reloaded == traj
        finally:
            monkey_patches.unpatch(paths)


    def test_export_trajectory_force(self, tmp_path):
        global paths
        outfile = tmp_path / "test_traj.db"
        outfile.touch()
        init_snap = make_1d_traj([0.0])[0]
        continue_condition = paths.LengthEnsemble(5).can_append
        traj = self.stupid.generate(init_snap, [continue_condition])
        assert len(traj) == 5
        assert outfile.exists()
        from openpathsampling.experimental.storage import (
            monkey_patches, Storage
        )
        paths = monkey_patches.monkey_patch_all(paths)

        try:
            self.stupid.export_trajectory(traj, outfile, force=True)
            assert outfile.exists()
            storage = Storage(outfile, mode='r')
            assert len(storage.trajectories) == 1
            reloaded = storage.trajectories[0]
            assert len(reloaded) == 5
            assert reloaded == traj
        finally:
            monkey_patches.unpatch(paths)

    def test_export_trajecory_file_exists_fail(self, tmp_path):
        global paths
        outfile = tmp_path / "test_traj.db"
        outfile.touch()
        init_snap = make_1d_traj([0.0])[0]
        continue_condition = paths.LengthEnsemble(5).can_append
        traj = self.stupid.generate(init_snap, [continue_condition])
        assert len(traj) == 5
        assert outfile.exists()
        from openpathsampling.experimental.storage import (
            monkey_patches, Storage
        )
        paths = monkey_patches.monkey_patch_all(paths)

        try:
            with pytest.raises(FileExistsError):
                self.stupid.export_trajectory(traj, outfile)
        finally:
            monkey_patches.unpatch(paths)

    def test_export_trajectory_custom_writer(self, tmp_path):
        # here we use a custom writer that will write to a string
        from openpathsampling.exports.trajectories.core import (
            TrajectoryWriter
        )
        outfile = tmp_path / "test_traj.txt"
        class StringXWriter(TrajectoryWriter):
            def _write(self, trajectory, filename):
                outstring = " ".join([
                    f"{snap.xyz[0][0]:.2f}" for snap in trajectory
                ])
                with open(filename, 'w') as f:
                    f.write(outstring)

        init_snap = make_1d_traj([0.0])[0]
        continue_condition = paths.LengthEnsemble(5).can_append
        traj = self.stupid.generate(init_snap, [continue_condition])
        assert len(traj) == 5
        writer = StringXWriter()
        self.stupid.export_trajectory(traj, outfile, writer=writer)
        with open(outfile, 'r') as f:
            outstring = f.read()

        assert outstring == ("0.00 " * 5)[:-1]
