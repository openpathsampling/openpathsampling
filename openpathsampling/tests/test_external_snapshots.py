import pytest

import openpathsampling as paths
import numpy as np

from openpathsampling.engines.external_snapshots.snapshot import (
    ExternalMDSnapshot, InternalizedMDSnapshot
)
from openpathsampling.engines.snapshot import SnapshotDescriptor

from openpathsampling.engines.external_engine import \
        _InternalizedEngineProxy

class MockEngine(object):
    SnapshotClass = ExternalMDSnapshot
    InternalizedSnapshotClass = InternalizedMDSnapshot
    def __init__(self, sequences, sleep_ms=0):
        self.sequences = sequences
        self.sleep_ms = sleep_ms
        self.descriptor = SnapshotDescriptor.construct(
            snapshot_class=ExternalMDSnapshot,
            snapshot_dimensions={'n_spatial': 2, 'n_atoms': 1}
        )
        self.internalized_engine = _InternalizedEngineProxy(self)

    def read_frame_data(self, filename, position):
        return self.sequences[filename][position]


class ErrorMockEngine(MockEngine):
    """Mock engine used to create the IndexError in load_details"""
    def __init__(self, sequences, sleep_ms=0):
        super(ErrorMockEngine, self).__init__(sequences, sleep_ms)
        self._sequences = self.sequences
        self.sequences = {k: [] for k in self.sequences}
        self.accessed = False

    def read_frame_data(self, filename, position):
        if self.accessed:
            self.sequences = self._sequences

        self.accessed = True
        return super(ErrorMockEngine, self).read_frame_data(filename,
                                                            position)


class TestExternalMDSnapshot(object):
    def setup(self):
        self.box = np.array([[1.0, 0.0], [0.0, 1.0]])
        self.vel = np.array([[1.0, 0.0]])


        self.engine = MockEngine(
            sequences={
                'foo': [(np.array([[0.0, 0.0]]), self.vel, self.box),
                        (np.array([[0.1, 0.0]]), self.vel, self.box),
                        (np.array([[0.2, 0.0]]), self.vel, self.box)]
            },
            sleep_ms=0.0
        )
        self.snapshots = [
            ExternalMDSnapshot(file_name='foo',
                               file_position=i,
                               engine=self.engine)
            for i in range(3)
        ]

    def test_init(self):
        for (i, snap) in enumerate(self.snapshots):
            assert snap._reversed is None
            assert snap.file_name == "foo"
            assert snap.file_position == i
            assert snap.engine == self.engine
            assert snap.velocity_direction == 1

    def test_load_details(self):
        for (i, snap) in enumerate(self.snapshots):
            snap.load_details()
            expected_xyz = np.array([[0.1 * i, 0.0]])
            expected_vel = np.array([[1.0, 0.0]])
            np.testing.assert_array_equal(snap._xyz, expected_xyz)
            np.testing.assert_array_equal(snap._velocities, expected_vel)
            np.testing.assert_array_equal(snap._box_vectors, self.box)

    def test_load_details_indexerror(self):
        engine = ErrorMockEngine(self.engine.sequences)
        snap = ExternalMDSnapshot(file_name='foo',
                                  file_position=0,
                                  engine=engine)
        snap.load_details()
        np.testing.assert_array_equal(snap._xyz, np.array([[0.0, 0.0]]))
        np.testing.assert_array_equal(snap._velocities, self.vel)
        np.testing.assert_array_equal(snap._box_vectors, self.box)

    def test_load_details_recursionerror(self):
        bad_snap = ExternalMDSnapshot(file_name='foo',
                                      file_position=3,
                                      engine=self.engine)
        with pytest.raises(RuntimeError):
            bad_snap.load_details()

    def test_set_details(self):
        xyz = np.array([[1.0, 1.0]])
        vel = np.array([[2.0, 2.0]])
        snap = ExternalMDSnapshot(file_name='bar', file_position=0)
        snap.set_details(xyz=xyz, velocities=vel, box_vectors=self.box)
        np.testing.assert_array_equal(snap._xyz, xyz)
        np.testing.assert_array_equal(snap._velocities, vel)
        np.testing.assert_array_equal(snap._box_vectors, self.box)

    def test_set_details_exists(self):
        snap = self.snapshots[0]
        with pytest.raises(RuntimeError):
            snap.set_details(xyz=np.array([[1.0, 1.0]]),
                             velocities=self.vel,
                             box_vectors=self.box)

    def test_clear_cache(self):
        snap = self.snapshots[0]
        snap.load_details()
        assert snap._xyz is not None
        assert snap._velocities is not None
        assert snap._box_vectors is not None
        snap.clear_cache()
        assert snap._xyz is None
        assert snap._velocities is None
        assert snap._box_vectors is None

    def test_reversed(self):
        traj = paths.Trajectory(self.snapshots)
        traj_rev = traj.reversed
        for idx in range(len(traj)):
            snap = traj[idx]
            snap_rev = traj_rev[-idx-1]
            np.testing.assert_array_equal(snap.xyz, snap_rev.xyz)
            np.testing.assert_array_equal(snap.box_vectors,
                                          snap_rev.box_vectors)
            assert snap.velocity_direction == 1
            assert snap_rev.velocity_direction == -1
            np.testing.assert_array_equal(snap.velocities,
                                          -snap_rev.velocities)

            assert snap._reversed == snap_rev
            assert snap_rev._reversed == snap

    def test_internalize(self):
        snap = self.snapshots[0]
        internal = snap.internalize()
        np.testing.assert_array_equal(snap.xyz, internal.xyz)
        np.testing.assert_array_equal(snap.velocities, internal.velocities)
        np.testing.assert_array_equal(snap.box_vectors, internal.box_vectors)

        # the way to do it for a trajectory
        traj_i = paths.Trajectory([s.internalize() for s in self.snapshots])
        traj_e = paths.Trajectory(self.snapshots)
        np.testing.assert_array_equal(traj_i.xyz, traj_e.xyz)
