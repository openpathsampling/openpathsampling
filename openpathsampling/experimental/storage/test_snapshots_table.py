from .snapshots_table import *
from ..simstore.serialization_helpers import get_uuid
import pytest
from unittest import mock

from collections import namedtuple

import numpy as np
from openpathsampling.engines import toy as toys


def make_engine():
    pes = toys.LinearSlope([0, 0], 0)
    topology = toys.Topology(n_spatial=2, masses=[1.0, 1.0], pes=pes)
    integ = toys.LeapfrogVerletIntegrator(dt=0.1)
    options = {'integ': integ, 'n_frames_max': 1000, 'n_steps_per_frame': 1}
    engine = toys.Engine(options=options, topology=topology)
    return engine

# mocks
UUIDHolder = namedtuple("UUIDHolder", ['uuid'])

class MockBackend(object):
    def __init__(self, tables):
        self.obj_tables = tables
        self.tables = {}
        self.uuid_map = {}
        self.saved = []
        self.update()

    def update(self):
        self.tables = {
            table_name: [UUIDHolder(get_uuid(obj)) for obj in table_list]
            for table_name, table_list in self.obj_tables.items()
        }
        self.uuid_map = {get_uuid(obj): obj
                         for table in self.obj_tables.values()
                         for obj in table}

    @property
    def table_to_class(self):
        return {table: toys.Snapshot for table in self.tables}

    def table_len(self, table):
        return len(self.obj_tables[table])

    def table_iterator(self, table):
        return iter(self.tables[table])

    def table_get_item(self, table, item):
        return self.tables[table][item]

    def load(self, uuids):
        return [self.uuid_map[uuid] for uuid in uuids]

    def save(self, snapshot):
        self.saved.append(snapshot)


class TestSnapshotsTable(object):
    def setup(self):
        self.engine_1 = make_engine()
        self.engine_2 = make_engine()
        snap_1_0 = toys.Snapshot(coordinates=np.array([0.0, 0.0]),
                                 velocities=np.array([0.0, 0.0]),
                                 engine=self.engine_1)
        snap_1_1 = toys.Snapshot(coordinates=np.array([1.0, 0.0]),
                                 velocities=np.array([1.0, 0.0]),
                                 engine=self.engine_1)
        snap_2_0 = toys.Snapshot(coordinates=np.array([0.0, 0.0]),
                                 velocities=np.array([0.0, 0.0]),
                                 engine=self.engine_2)
        self.all_snaps = [snap_1_0, snap_1_1, snap_2_0]

        backend = MockBackend({'snapshot0': [snap_1_0, snap_1_1],
                               'snapshot1': [snap_2_0]})
        self.storage = mock.Mock(backend=backend, load=backend.load,
                                 save=backend.save)

        self.snapshots = SnapshotsTable(self.storage)

    def test_initialization(self):
        assert len(self.snapshots.tables) == 2

    def test_iter(self):
        for truth, beauty in zip(self.all_snaps, self.snapshots):
            assert truth == beauty

    def test_getitem(self):
        for num, obj in enumerate(self.all_snaps):
            assert self.snapshots[num] == obj

    def test_len(self):
        assert len(self.snapshots) == 3

    @pytest.mark.parametrize('engine_num', [1, 2])
    def test_snapshots_for_engine(self, engine_num):
        engine = {1: self.engine_1, 2: self.engine_2}[engine_num]
        snaps = {1: self.all_snaps[:2], 2: self.all_snaps[2:]}[engine_num]
        snapshot_table = self.snapshots.snapshots_for_engine(engine)
        for truth, beauty in zip(snaps, snapshot_table):
            assert truth == beauty

    def test_cache_all(self):
        pytest.skip()

    def test_save(self):
        new_snap = toys.Snapshot(coordinates=np.array([-1.0, -1.0]),
                                 velocities=np.array([0.0, 0.0]),
                                 engine=make_engine())
        assert self.storage.backend.saved == []  # mock
        self.snapshots.save(new_snap)
        assert self.storage.backend.saved == [new_snap]
