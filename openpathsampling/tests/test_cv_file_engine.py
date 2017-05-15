from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises, assert_true)

from nose.plugins.skip import Skip, SkipTest

from test_helpers import data_filename, assert_items_almost_equal

import openpathsampling as paths
import openpathsampling.engines.cv_file as cv_engine
import openpathsampling.engines.toy as toys
import os
import numpy as np

import logging
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

class testCVFileEngine(object):
    def setUp(self):
        self.engine = cv_engine.Engine(
            options={
                'cv_names': ['alpha', 'beta', 'gamma']
            }
        )
        self.filename = data_filename('cv_engine_trajectory.data')

    def tearDown(self):
        pass

    def test_initialization_as_spatial(self):
        engine = cv_engine.Engine(
            options={
                'cv_names': ['alpha', 'beta', 'gamma'],
                'n_spatial': 1
            }
        )
        assert_equal(engine.n_atoms, 3)
        assert_equal(engine.n_spatial, 1)

    def test_initialization_no_options(self):
        engine = cv_engine.Engine()
        assert_equal(engine.n_atoms, 1)
        assert_equal(engine.n_spatial, 1)
        assert_equal(len(engine.cv), 1)
        assert_equal(engine.cv.keys(), ["x"])
        assert_equal(engine.cv['x'].name, "x")

    def test_trajectory_from_file(self):
        trajectory = self.engine.trajectory_from_file(self.filename)
        assert_items_almost_equal(self.engine.cv['alpha'](trajectory),
                                  [1.0, 3.0, 4.0, 8.0])
        assert_items_almost_equal(self.engine.cv['beta'](trajectory),
                                  [2.0, 3.0, 5.0, 9.0])
        assert_items_almost_equal(self.engine.cv['gamma'](trajectory),
                                  [3.0, 3.0, 6.0, 0.0])

    def test_read_frame_from_file(self):
        frame0 = self.engine.read_frame_from_file(self.filename, 0)
        frame1 = self.engine.read_frame_from_file(self.filename, 1)
        assert_equal(np.all(frame0.coordinates == np.array([1.0, 2.0, 3.0])),
                     True)
        assert_equal(np.all(frame1.coordinates == np.array([3.0, 3.0, 3.0])),
                     True)

    def test_write_frame_to_file(self):
        coords = [6.0, 8.0, 9.0]
        snap = toys.Snapshot(coordinates=np.array([coords]),
                             velocities=np.array([[0.0, 0.0, 0.0]]))
        self.engine.write_frame_to_file(data_filename("frame.out"), snap,
                                        mode="w")
        f = open(data_filename("frame.out"))
        lines = [l for l in f]
        assert_equal(len(lines), 1)
        splitted = map(float, lines[0].split())
        assert_items_almost_equal(splitted, coords)
        if os.path.isfile(data_filename("frame.out")):
            os.remove(data_filename("frame.out"))
