from builtins import object
import openpathsampling.engines.toy as toys
import numpy as np

from nose.tools import assert_equal, assert_almost_equal
from nose.plugins.skip import SkipTest
import logging

logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

class TestToySnapshotFeatures(object):
    def setup(self):
        integ = toys.LangevinBAOABIntegrator(dt=0.002,
                                             temperature=1.0,
                                             gamma=2.5)
        options={'integ': integ, 'n_frames_max': 5}
        topology_2D = toys.Topology(n_spatial=2,
                                    masses=[1.0],
                                    pes=None)
        engine_2D = toys.Engine(options=options, topology=topology_2D)
        self.snap_2D = toys.Snapshot(coordinates=np.array([[0.0, 0.0]]),
                                     velocities=np.array([[2.0, 4.0]]),
                                     engine=engine_2D)

        topology_6D = toys.Topology(n_spatial=3,
                                    n_atoms=2,
                                    masses=[1.0, 2.0],
                                    pes=None)
        engine_6D = toys.Engine(options=options, topology=topology_6D)
        self.snap_6D = toys.Snapshot(coordinates=np.array([[0.0, 0.0, 0.0],
                                                           [0.0, 0.0, 0.0]]),
                                     velocities=np.array([[1.0, 0.0, 0.0],
                                                          [1.0, 0.0, 0.0]]),
                                     engine=engine_6D)

    def test_n_degrees_of_freedom(self):
        assert_equal(self.snap_2D.n_degrees_of_freedom, 2)
        assert_equal(self.snap_6D.n_degrees_of_freedom, 6)

    def test_instantaneous_temperature(self):
        # KE = 0.5 * (1.0*2.0**2 + 1.0*4.0**2) = 10.0
        # T = 2 * KE / 2 = KE = 10.0
        assert_almost_equal(self.snap_2D.instantaneous_temperature, 10.0)
        # KE = 0.5 * (1.0*1.0**2 + 2.0*1.0**2) = 1.5
        # T = 2 * KE / 6 = 0.5
        assert_almost_equal(self.snap_6D.instantaneous_temperature, 0.5)
